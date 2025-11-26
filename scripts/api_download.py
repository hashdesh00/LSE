# scripts/api_download.py

from __future__ import annotations

import logging
from pathlib import Path
from typing import Dict, Optional, List, Any

import pandas as pd
import numpy as np
import requests
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry
from datetime import datetime

import pyarrow as pa
import pyarrow.parquet as pq

# -----------------------------------------------------------------------------------------------------------------------------
# Configuration
# -----------------------------------------------------------------------------------------------------------------------------

# Base paths
PROJECT_ROOT = Path(__file__).resolve().parent.parent
DATA_RAW_DIR = PROJECT_ROOT / "data_raw"
PORTWATCH_RAW_DIR = DATA_RAW_DIR / "portwatch"
OFFICIAL_RAW_DIR = DATA_RAW_DIR /"official"

PORTWATCH_DAILY_TRADE_URL = ("https://services9.arcgis.com/weJ1QsnbMYJlCHdG/arcgis/rest/services/""Daily_Trade_Data/FeatureServer/0/query")
PORTWATCH_PORTS_DB_URL = ("https://services9.arcgis.com/weJ1QsnbMYJlCHdG/ArcGIS/rest/services/PortWatch_ports_database/FeatureServer/query")
PORTWATCH_SPILLOVER_QUERY = ("https://services9.arcgis.com/weJ1QsnbMYJlCHdG/ArcGIS/rest/services/spillovers_country/FeatureServer/query")

PORTWATCH_RAW_DIR.mkdir(parents=True, exist_ok=True)
OFFICIAL_RAW_DIR.mkdir(parents=True, exist_ok=True)


# If downloading official trade via HTTP, set the URL here
FOREIGN_TRADE_TIMESERIES_URL : Optional[str] = None          # fill if needed

# Logging setup
logger = logging.getLogger(__name__)
if not logger.handlers:
    handler = logging.StreamHandler()
    formatter = logging.Formatter("[%(asctime)s] %(levelname)s - %(name)s - %(message)s",
                                  datefmt = "%Y-%m-%d %H:%M:%S"
    )
    handler.setFormatter(formatter)
    logger.addHandler(handler)
    logger.setLevel(logging.INFO)


# ---------------------------------------------------------------------------
# HTTP session with retries
# ---------------------------------------------------------------------------

def create_session(
    total_retries: int = 5,
    backoff_factor: float = 0.5,
    status_forcelist: Optional[List[int]] = None,
) -> requests.Session:
    """
    Create a shared requests.Session with retry logic.

    Parameters
    ----------
    total_retries : int
        Total number of retries for failed requests.
    backoff_factor : float
        Sleep factor between retries (exponential backoff).
    status_forcelist : list[int], optional
        HTTP status codes that should trigger a retry.

    Returns
    -------
    requests.Session
    """
    if status_forcelist is None:
        status_forcelist = [429, 500, 502, 503, 504]

    retry = Retry(
        total=total_retries,
        read=total_retries,
        connect=total_retries,
        backoff_factor=backoff_factor,
        status_forcelist=status_forcelist,
        allowed_methods=["GET", "POST"],
        raise_on_status=False,
    )

    adapter = HTTPAdapter(max_retries=retry)

    session = requests.Session()
    session.mount("http://", adapter)
    session.mount("https://", adapter)

    return session


# Create one global session for all API calls in this module
SESSION = create_session()


# ---------------------------------------------------------------------------
# Core helper for ArcGIS FeatureServer
# ---------------------------------------------------------------------------

def fetch_arcgis_features(
    url: str,
    where: str,
    out_fields: str = "*",
    result_record_count: int = 2000,
    session: Optional[requests.Session] = None,
    extra_params: Optional[Dict[str, Any]] = None,
) -> pd.DataFrame:
    """
    Fetch all features from an ArcGIS FeatureServer layer using pagination.

    Parameters
    ----------
    url : str
        ArcGIS FeatureServer 'query' endpoint.
    where : str
        WHERE clause (e.g. "country = 'Malaysia' AND year >= 2019").
    out_fields : str
        Fields to return, default "*".
    result_record_count : int
        Page size for pagination.
    session : requests.Session, optional
        HTTP session to use (with retries).
    extra_params : dict, optional
        Additional parameters to pass to the query.

    Returns
    -------
    pd.DataFrame
        DataFrame of all features.
    """
    if session is None:
        session = SESSION

    params = {
        "where": where,
        "outFields": out_fields,
        "outSR": "4326",
        "f": "json",
        "resultRecordCount": result_record_count,
        "resultOffset": 0,
        "returnExceededLimitFeatures": True,
    }

    if extra_params:
        params.update(extra_params)

    all_records: List[Dict[str, Any]] = []
    logger.info(f"Fetching ArcGIS data from: {url}")
    logger.info(f"WHERE clause: {where}")

    while True:
        logger.debug(f"Requesting offset {params['resultOffset']}")
        resp = session.get(url, params=params, timeout=60)
        try:
            resp.raise_for_status()
        except requests.HTTPError as e:
            logger.error(f"HTTP error: {e}")
            logger.error(f"Response text: {resp.text[:500]}")
            raise

        data = resp.json()

        if "error" in data:
            logger.error(f"ArcGIS error: {data['error']}")
            raise RuntimeError(f"ArcGIS error: {data['error']}")

        features = data.get("features", [])
        if not features:
            logger.info("No more features returned from server.")
            break

        # features are list of {"attributes": {...}, "geometry": {...}}
        for feat in features:
            attrs = feat.get("attributes", {})
            all_records.append(attrs)

        logger.info(f"Fetched {len(features)} records (total so far: {len(all_records)})")

        if not data.get("exceededTransferLimit", False):
            # No more pages, server confirms
            logger.info("Server indicates no more features are available (exceededTransferLimit=false).")
            break

        # Increment offset for next page
        params["resultOffset"] += len(features)

    if not all_records:
        logger.warning("No records fetched; returning empty DataFrame.")
        return pd.DataFrame()

    df = pd.DataFrame(all_records)
    logger.info(f"Final DataFrame shape: {df.shape}")
    return df

# ---------------------------------------------------------------------------
# PortWatch-specific extractors
# ---------------------------------------------------------------------------

def build_year_where_clause(
    base_filters: Optional[List[str]] = None,
    start_year: Optional[int] = None,
    end_year: Optional[int] = None,
) -> str:
    """
    Helper to build a WHERE clause with optional year filters.
    """
    clauses: List[str] = base_filters.copy() if base_filters else []

    if start_year is not None:
        clauses.append(f"year >= {start_year}")
    if end_year is not None:
        clauses.append(f"year <= {end_year}")

    if not clauses:
        return "1=1"  # ArcGIS 'no filter'
    return " AND ".join(clauses)



# ---------------------------------------------------------------------------
# Official trade downloader (optional HTTP)
# ---------------------------------------------------------------------------

def download_official_foreign_trade(
    force: bool = False,
    session: Optional[requests.Session] = None,
) -> Optional[Path]:
    """
    (Optional) Download official foreign trade timeseries via HTTP and cache it.

    If your current approach uses a SharePoint helper inside the Bank network
    (e.g. read_sharepoint_file), you might not need this HTTP version. This is
    here as a template.

    Returns
    -------
    Path or None
        Path to saved parquet if URL is configured, else None.
    """
    if FOREIGN_TRADE_TIMESERIES_URL is None:
        logger.warning("FOREIGN_TRADE_TIMESERIES_URL is not set; skipping download.")
        return None

    session = session or SESSION
    out_path = OFFICIAL_RAW_DIR / "foreign_trade_timeseries.parquet"

    if out_path.exists() and not force:
        logger.info(f"Using cached official trade file: {out_path}")
        return out_path

    logger.info(f"Downloading official trade data from: {FOREIGN_TRADE_TIMESERIES_URL}")
    resp = session.get(FOREIGN_TRADE_TIMESERIES_URL, timeout=60)
    try:
        resp.raise_for_status()
    except requests.HTTPError as e:
        logger.error(f"HTTP error when downloading official trade: {e}")
        logger.error(f"Response text: {resp.text[:500]}")
        raise

    # Example: if server returns CSV; adjust as needed
    df = pd.read_csv(pd.compat.StringIO(resp.text))  # or pd.read_parquet if binary
    logger.info(f"Saving official trade data to {out_path}")
    df.to_parquet(out_path, index=False)
    return out_path

# ---------------------------------------------------------------------------
# Daily trade download (year-chunked)
# ---------------------------------------------------------------------------

def download_portwatch_daily_trade(
    country: str,
    start_year: int = 2019,
    end_year: int | None = None,
    force: bool = False,
) -> Path:
    """
    Download Daily_Trade_Data for a given country across multiple years.

    Saves:
        data_raw/portwatch/portwatch_daily_trade_{country}_from{start}_to{end}.parquet
    """
    if end_year is None:
        end_year = datetime.utcnow().year

    safe_country = country.replace(" ", "_")
    filename = f"portwatch_daily_trade_{safe_country}_from{start_year}_to{end_year}.parquet"
    out_path = PORTWATCH_RAW_DIR / filename

    if out_path.exists() and not force:
        logger.info(f"[Daily_Trade] Using cached file: {out_path}")
        return out_path

    logger.info(f"[Daily_Trade] Downloading PortWatch for {country} ({start_year}–{end_year})")

    all_chunks: list[pd.DataFrame] = []

    for year in range(start_year, end_year + 1):
        where = f"country = '{country}' AND year = {year}"
        logger.info(f"  • Querying year {year} with WHERE = {where}")

        df_year = fetch_arcgis_features(
            url=PORTWATCH_DAILY_TRADE_URL,
            where=where,
            out_fields="*",
            result_record_count=2000,
            session=SESSION,        # <-- USE EXISTING GLOBAL SESSION
        )

        if not df_year.empty:
            all_chunks.append(df_year)
        else:
            logger.warning(f"  • No data returned for year {year}")

    if not all_chunks:
        raise RuntimeError(f"No Daily_Trade_Data returned for {country} {start_year}–{end_year}")

    df = pd.concat(all_chunks, ignore_index=True)
    logger.info(f"[Daily_Trade] Final rows: {df.shape[0]:,}")

    df.to_parquet(out_path, index=False)
    return out_path


# ---------------------------------------------------------------------------
# Ports database download (metadata, shares, vessel counts)
# ---------------------------------------------------------------------------

def download_portwatch_ports_metadata(
    country: str,
    force: bool = False,
) -> Path:
    """
    Download PortWatch_ports metadata for a country.

    Output: data_raw/portwatch/portwatch_ports_metadata_{country}.parquet
    """
    safe_country = country.replace(" ", "_")
    filename = f"portwatch_ports_metadata_{safe_country}.parquet"
    out_path = PORTWATCH_RAW_DIR / filename

    if out_path.exists() and not force:
        logger.info(f"[Ports metadata] Using cached file: {out_path}")
        return out_path

    where = f"country = '{country}'"
    logger.info(f"[Ports metadata] Downloading with WHERE = {where}")

    df = fetch_arcgis_features(
        url=PORTWATCH_PORTS_DB_URL,
        where=where,
        out_fields="*",
        result_record_count=2000,
    )
    if df.empty:
        raise RuntimeError(f"No ports metadata returned for {country}")

    df.to_parquet(out_path, index=False)
    return out_path
