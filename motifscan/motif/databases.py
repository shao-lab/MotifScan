"""
motifscan.motif.databases
-------------------------

Get motifs from public databases.
"""

import logging
import os
import re
import sys
from datetime import datetime

import requests

from motifscan.exceptions import RemoteMotifPFMsNotFoundError

logger = logging.getLogger(__name__)


class JasparDatabase:
    """JASPAR motif database for the JASPAR latest 2020 release, including
    JASPAR CORE (non-redundant/redundant set) and JASPAR collections.

    Attributes
    ----------
    name : str
        Database name, `JASPAR2020` here.
    core_taxons : list of str
        Taxons of JASPAR CORE collection.
    other_collections : list of str
        JASPAR other collections.
    url_core_fmt : str
        URL formatter for motif sets in JASPAR CORE.
    url_other_collections_fmt : str
        URL formatter for motif sets in JASPAR other collections.

    Notes
    -----
        - JASPAR CORE non-redundant set should be used in most situations.
        - JASPAR CORE redundant set includes PFMs released in earlier versions.
    """

    def __init__(self):
        self.name = 'JASPAR2020'
        self.core_taxons = ['vertebrates', 'plants', 'insects', 'nematodes',
                            'fungi', 'urochordates']
        self.other_collections = ['CNE', 'PHYLOFACTS', 'SPLICE', 'POLII',
                                  'FAM', 'PBM', 'PBM_HOMEO', 'PBM_HLH',
                                  'UNVALIDATED']
        url_prefix = "http://jaspar.genereg.net/download/"
        self.url_core_fmt = \
            url_prefix + "CORE/JASPAR2020_CORE_{0}_pfms_jaspar.txt"
        self.url_other_collections_fmt = \
            url_prefix + "collections/JASPAR2020_{0}_pfms_jaspar.txt"
        self._pfms_core = None
        self._pfms_other_collections = self.other_collections

    @property
    def pfms_core(self):
        """Returns available motif PFMs sets in JASPAR CORE.

        Returns
        -------
        list of str
            List of available motif PFMs sets in JASPAR CORE.
        """
        if self._pfms_core is None:
            pfms = []
            for taxon in self.core_taxons:
                pfms.append(f"{taxon}_non-redundant")
                pfms.append(f"{taxon}_redundant")
            self._pfms_core = pfms
        return self._pfms_core

    @property
    def pfms_other_collections(self):
        """Returns available motif PFMs sets in JASPAR other collections.

        Returns
        -------
        list of str
            List of available motif PFMs sets in JASPAR other collections.
        """
        return self._pfms_other_collections

    @staticmethod
    def _download_pfms(pfms_url, download_dir):
        """Download PFMs from `pfms_url` and save to `download_dir`."""
        if not os.path.isdir(os.path.dirname(download_dir)):
            os.makedirs(download_dir)
        base_name = os.path.basename(pfms_url)
        dst = os.path.join(download_dir, base_name)
        try:
            logger.debug(f"Downloading {pfms_url}")
            r = requests.get(pfms_url, stream=True)
            r.raise_for_status()
            with open(dst, 'wb') as f:
                for chunk in r.iter_content(chunk_size=1024):
                    f.write(chunk)
        except requests.HTTPError as e:
            logger.error(f"Failed to download due to an HTTPError: {e}")
            sys.exit(1)
        return dst

    @staticmethod
    def _write_readme(database, pfms_name, download_dir):
        logger.debug("Writing the README file")
        readme_file = os.path.join(download_dir, "README")
        time_now = datetime.now().strftime("%Y-%m-%d %H:%M")
        with open(readme_file, 'w') as f_out:
            f_out.write(f"{database}\t{pfms_name}\tDownloaded at {time_now}\n")

    def download_core(self, pfms_name, download_dir):
        """Download a motif PFMs set from JASPAR CORE.

        Parameters
        ----------
        pfms_name : str
            Motif PFMs set name.
        download_dir : str
            Directory to write downloaded file.

        Returns
        -------
        dst_pfms : str
            Path of saved PFMs file

        Raises
        ------
        RemotePfmsNotFoundError
            If the specified PFMs is not in JASPAR CORE.
        """
        m = re.match(r"^([a-z]+)_(non-)?redundant$", pfms_name)
        if m and m.group(1) in self.core_taxons:
            logger.info(
                f"Downloading motif PFMs set {pfms_name!r} from JASPAR CORE")
            pfms_url = self.url_core_fmt.format(pfms_name)
            dst_pfms = self._download_pfms(pfms_url, download_dir)
            self._write_readme("JASPAR2020_CORE", pfms_name, download_dir)
        else:
            raise RemoteMotifPFMsNotFoundError("JASPAR CORE", pfms_name)
        return dst_pfms

    def download_other_collections(self, pfms_name, download_dir):
        """Download motif PFMs from JASPAR collections (other than CORE).

        Parameters
        ----------
        pfms_name : str
            Motif PFMs name.
        download_dir : str
            Directory to write downloaded file.

        Returns
        -------
        dst_pfms : str
            Path of saved PFMs file

        Raises
        ------
        RemotePfmsNotFoundError
            If the specified PFMs is not in JASPAR collections.
        """
        m = re.match(r"^([A-Z_]+)", pfms_name)
        if m and m.group(1) in self.other_collections:
            logger.info(f"Downloading motif PFMs set {pfms_name!r} from "
                        f"JASPAR Collections")
            pfms_url = self.url_other_collections_fmt.format(pfms_name)
            dst_pfms = self._download_pfms(pfms_url, download_dir)
            self._write_readme("JASPAR2020_Collections", pfms_name,
                               download_dir)
        else:
            raise RemoteMotifPFMsNotFoundError("JASPAR Collections", pfms_name)
        return dst_pfms

    @staticmethod
    def get_motif_info(matrix_id):
        """Get detailed information about a JASPAR motif.

        Parameters
        ----------
        matrix_id : str
            The matrix id of the JASPAR motif.

        Returns
        -------
        motif_info : dict
            A dict object contains the motif information.
        """
        motif_info = {}
        try:
            logger.debug(f"Getting motif info: {matrix_id}")
            url = f"http://jaspar.genereg.net/api/v1/matrix/{matrix_id}/"
            r = requests.get(url)
            r.raise_for_status()
            motif_info = r.json()
        except requests.HTTPError as e:
            logger.error(f"Failed to get motif info from JASPAR: {e}")
        return motif_info
