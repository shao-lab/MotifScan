"""
motifscan.io.utils
------------------

Utilities for file operations.
"""

import re
import gzip
import shutil
import tarfile
import zipfile


def replace_special_char(name):
    return re.sub('[-:./*]', '_', name)


def extract_zip(src, dst):
    """Extract and merge zip compressed (.zip) file.

    Parameters
    ----------
    src : str
        Path of zip compressed file.
    dst : str
        Path to write the decompressed file.
    """
    with zipfile.ZipFile(src, 'r') as fin, open(dst, 'wb') as fout:
        for member in fin.infolist():
            shutil.copyfileobj(fin.open(member), fout)


def extract_gzip(src, dst):
    """Extract gzip compressed (.gz) file.

    Parameters
    ----------
    src : str
        Path of gzip compressed file.
    dst : str
        Path to write the decompressed file.
    """
    with gzip.open(src, 'rb') as fin, open(dst, 'wb') as fout:
        shutil.copyfileobj(fin, fout)


def extract_targz(src, dst):
    """Extract and merge gzip compressed tar archive (.tar.gz) file.

    Parameters
    ----------
    src : str
        Path of gzip compressed tar archive file.
    dst : str
        Path to write the decompressed file.
    """
    with tarfile.open(src, 'r:gz') as fin, open(dst, 'wb') as fout:
        for member in fin.getmembers():
            if member.isfile():
                shutil.copyfileobj(fin.extractfile(member), fout)


def copy_file(src, dst):
    """Copy file.

    Parameters
    ----------
    src : str
        Path of the file to copy from.
    dst : str
        Path to write the copied file.
    """
    with open(src, 'rb') as fin, open(dst, 'wb') as fout:
        shutil.copyfileobj(fin, fout)


def merge_files(src, dst):
    """Merge(concatenate) files into one single file.

    Parameters
    ----------
    src : str
        Paths of the files to merge from.
    dst : str
        Path to write the merged file.
    """
    with open(dst, 'wb') as fout:
        for tmp_file in src:
            with open(tmp_file, 'rb') as fin:
                shutil.copyfileobj(fin, fout)


def merge_extracted_files(src, dst):
    """Extract compressed files and merge them into a single file.

    Parameters
    ----------
    src : str
        Path of compressed files.
    dst : str
        Path to write the merged file.
    """
    if src.endswith('tar.gz'):
        extract_targz(src, dst)
    elif src.endswith('zip'):
        extract_zip(src, dst)
    elif src.endswith('gz'):
        extract_gzip(src, dst)
    else:
        raise TypeError(f"unsupported compressed file: {src}")
