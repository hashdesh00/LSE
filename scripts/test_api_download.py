
import unittest
from unittest.mock import patch, MagicMock
import pandas as pd
from scripts.api_download import fetch_arcgis_features

class TestApiDownload(unittest.TestCase):

    @patch('scripts.api_download.requests.Session')
    def test_fetch_arcgis_features_pagination(self, mock_session):
        # Mock the API response
        mock_response_page1 = MagicMock()
        mock_response_page1.json.return_value = {
            "features": [{"attributes": {"id": 1}}, {"attributes": {"id": 2}}],
            "exceededTransferLimit": True
        }
        mock_response_page1.raise_for_status.return_value = None

        mock_response_page2 = MagicMock()
        mock_response_page2.json.return_value = {
            "features": [{"attributes": {"id": 3}}],
            "exceededTransferLimit": False
        }
        mock_response_page2.raise_for_status.return_value = None

        mock_session.get.side_effect = [mock_response_page1, mock_response_page2]

        # Call the function
        df = fetch_arcgis_features("http://fake-url.com", "1=1", session=mock_session)

        # Assertions
        self.assertEqual(len(df), 3)
        self.assertEqual(mock_session.get.call_count, 2)
        self.assertEqual(df['id'].tolist(), [1, 2, 3])

if __name__ == '__main__':
    unittest.main()
