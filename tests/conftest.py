import os
import shutil

import pytest

from motifscan.config import Config


@pytest.fixture(scope='session')
def data_dir():
    """Return the directory of data files."""
    return os.path.join(os.path.dirname(os.path.abspath(__file__)), 'data')


@pytest.fixture(scope='session')
def tmp_dir(data_dir):
    """Generate a temporary output directory to write output files."""
    tmp_dir = os.path.join(data_dir, 'motifscan_tmp_output')
    if os.path.isdir(tmp_dir):
        shutil.rmtree(tmp_dir)
    os.makedirs(tmp_dir)
    yield tmp_dir
    if os.path.isdir(tmp_dir):
        shutil.rmtree(tmp_dir)


@pytest.fixture(scope='function')
def config(tmp_dir):
    """Returns a config instance to test the configuration of MotifScan."""
    return Config(os.path.join(tmp_dir, '.motifscanrc'))


@pytest.fixture(scope='session')
def genome_root(data_dir):
    """Returns the path of test genome."""
    return os.path.join(data_dir, 'genomes')


@pytest.fixture(scope='session')
def motif_root(data_dir):
    """Returns the path of test motif set."""
    return os.path.join(data_dir, 'motifs')


@pytest.fixture(scope='session')
def region_root(data_dir):
    """Returns the path of test regions."""
    return os.path.join(data_dir, 'regions')
