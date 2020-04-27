"""
motifscan.genome.databases
--------------------------

Get genome assemblies from public databases.
"""

import logging
import os
import sys
from collections import namedtuple
from xml.etree import ElementTree

import requests
from tqdm import tqdm

from motifscan.exceptions import RemoteGenomeNotFoundError, \
    RemoteGenomeFileNotFoundError

logger = logging.getLogger(__name__)

Assembly = namedtuple('Assembly', ['id', 'description'])


class UcscDatabase:
    """UCSC genome database.

    Attributes
    ----------
    name : str
        Database name, `UCSC` here.
    info_page : str
        The DAS info URL of the UCSC database.
    seq_url_fmts : list of str
        URL patterns for the sequence file.
    gene_url_fmts : list of str
        URL patterns for the gene annotation file.
    """

    def __init__(self):
        self.name = 'UCSC'
        self._assemblies = None
        self.info_page = "https://genome.ucsc.edu/cgi-bin/das/dsn"
        url_prefix = "http://hgdownload.soe.ucsc.edu/goldenPath/"
        self.seq_url_fmts = [
            url_prefix + "{0}/bigZips/{1}.fa.gz",
            url_prefix + "{0}/bigZips/chromFa.tar.gz",
            url_prefix + "{0}/bigZips/{1}.chromFa.tar.gz",
            url_prefix + "{0}/bigZips/chromFa.zip",
        ]
        self.gene_url_fmts = [url_prefix + "{0}/database/refGene.txt.gz"]

    @property
    def assemblies(self):
        """Returns available genome assemblies in the UCSC database.

        Returns
        -------
        list of `Assembly`
            Available assemblies in the UCSC database. Items are namedtuple
            objects with `id` and `description` as attributes.
        """
        if self._assemblies is None:
            r = requests.get(self.info_page)
            root = ElementTree.fromstring(r.text)
            assemblies = []
            for child in root:
                if child[0].tag == 'SOURCE' and child[1].tag == 'DESCRIPTION':
                    assemblies.append(Assembly(id=child[0].attrib['id'],
                                               description=child[1].text))
            self._assemblies = assemblies
        return self._assemblies

    def search(self, keyword):
        """Search genome assemblies by matching the given keyword either in
        the id field or the description field.

        Parameters
        ----------
        keyword : str
            The keyword to search for (case-insensitive).

        Yields
        ------
        assembly : `Assembly`
            Assemblies matching the specified `keyword`.
        """
        for assembly in self.assemblies:
            if (keyword.lower() in assembly.id.lower() or
                    keyword.lower() in assembly.description.lower()):
                yield assembly

    def get_sequence_url(self, assembly):
        """Returns the sequence file URL of the specified genome assembly.

        Parameters
        ----------
        assembly : str
            Genome assembly name (case-sensitive).

        Returns
        -------
        url : str
            The URL of the sequence file.

        Raises
        ------
        RemoteGenomeNotFoundError
            If the specified assembly is not in the database.
        RemoteGenomeFileNotFoundError
            If the sequence file does not exists or is not accessible.
        """
        found_assembly = False
        for tmp_assembly in self.assemblies:
            if tmp_assembly.id == assembly:
                found_assembly = True
                for pattern in self.seq_url_fmts:
                    url = pattern.format(assembly, assembly)
                    logger.debug(f"Trying URL: {url}")
                    r = requests.head(url)
                    if r.status_code == requests.codes.ok:
                        logger.debug("Succeed")
                        return url
                    else:
                        logger.debug(
                            f"Failed with status code: {r.status_code}")
        if not found_assembly:
            raise RemoteGenomeNotFoundError(self.name, assembly)
        else:
            raise RemoteGenomeFileNotFoundError(self.name, assembly,
                                                'sequence')

    def get_gene_url(self, assembly):
        """Returns the gene annotation file URL of the specified genome
        assembly.

        Parameters
        ----------
        assembly : str
            Genome assembly name (case-sensitive).

        Returns
        -------
        url : str
            The URL of the gene annotation file.

        Raises
        ------
        RemoteGenomeNotFoundError
            If the specified assembly is not in the database.
        RemoteGenomeFileNotFoundError
            If the gene annotation file does not exists or is not accessible.
        """
        found_assembly = False
        for tmp_assembly in self.assemblies:
            if tmp_assembly.id == assembly:
                found_assembly = True
                for pattern in self.gene_url_fmts:
                    url = pattern.format(assembly)
                    logger.debug(f"Trying url: {url}")
                    r = requests.head(url)
                    if r.status_code == requests.codes.ok:
                        logger.debug("Succeed")
                        return url
                    else:
                        logger.debug(
                            f"Failed with status code: {r.status_code}")
        if not found_assembly:
            raise RemoteGenomeNotFoundError(self.name, assembly)
        else:
            raise RemoteGenomeFileNotFoundError(self.name, assembly,
                                                'annotation')

    @staticmethod
    def _download_ucsc_file(url, download_dir):
        """Download file from `url` and save into `download_dir`."""
        if not os.path.isdir(download_dir):
            os.makedirs(download_dir)
        base_name = os.path.basename(url)
        dst = os.path.join(download_dir, base_name)
        try:
            logger.debug(f"Downloading {url}")
            r = requests.get(url, stream=True)
            r.raise_for_status()
            total_size = int(r.headers.get('Content-Length'))
            with open(dst, 'wb') as f:
                with tqdm(total=total_size, unit='B', unit_scale=True,
                          desc=base_name) as pbar:
                    for chunk in r.iter_content(chunk_size=1024):
                        f.write(chunk)
                        pbar.update(len(chunk))
        except requests.HTTPError as e:
            logger.error(f"Failed to download due to an HTTPError: {e}")
            sys.exit(1)
        return dst

    def download_sequence(self, assembly, download_dir):
        """Download the sequence file of the specified genome assembly.

        Parameters
        ----------
        assembly : str
            Genome assembly name (case-sensitive).
        download_dir : str
            Directory to write the downloaded file.

        Returns
        -------
        dst : str
            Path of saved sequence file.
        """
        url = self.get_sequence_url(assembly)
        logger.info("Downloading the sequence file")
        dst = self._download_ucsc_file(url, download_dir)
        return dst

    def download_gene(self, assembly, download_dir):
        """Download the gene annotation file of the specified genome assembly.

        Parameters
        ----------
        assembly : str
            Genome assembly name (case-sensitive).
        download_dir : str
            Directory to write the downloaded file.

        Returns
        -------
        dst : str
            Path of saved gene annotation file.
        """
        url = self.get_gene_url(assembly)
        logger.info("Downloading the gene annotation file")
        dst = self._download_ucsc_file(url, download_dir)
        return dst
