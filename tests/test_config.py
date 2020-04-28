import os

import pytest

from motifscan.config import Config, user_rc_path, user_genome_dir, \
    user_motif_dir
from motifscan.exceptions import InvalidConfigFileError, GenomeNotFoundError, \
    MotifSetNotFoundError


def test_invalid_config(data_dir):
    with pytest.raises(InvalidConfigFileError):
        Config(os.path.join(data_dir, 'invalid.motifscanrc'))


def test_config_init(config):
    assert len(config._config.sections()) == 3
    assert config._config.has_section('motifscan')
    assert config._config.has_section('genome')
    assert config._config.has_section('motif')
    assert config._config.get('motifscan', 'genome_dir') == user_genome_dir
    assert config._config.get('motifscan', 'motif_dir') == user_motif_dir
    config = Config(path=None)
    assert config.path == user_rc_path


def test_get_genome_dir(config):
    assert config.get_genome_dir() == user_genome_dir


def test_set_genome_dir(config):
    path = 'test_dir'
    config.set_genome_dir(path)
    assert config.get_genome_dir() == path


def test_get_motif_dir(config):
    assert config.get_motif_dir() == user_motif_dir


def test_set_motif_dir(config):
    path = 'test_dir'
    config.set_motif_dir(path)
    assert config.get_motif_dir() == path


def test_set_genome_path(config):
    name = 'hg19'
    path = 'test_dir'
    config.set_genome_path(name, path)
    assert config._config.get('genome', name) == path


def test_get_genome_path(config):
    name = 'hg19'
    path = 'test_dir'
    config.set_genome_path(name, path)
    assert config.get_genome_path(name) == path
    with pytest.raises(GenomeNotFoundError):
        assert config.get_genome_path('mm9')


def test_remove_genome_path(config):
    name = 'hg19'
    path = 'test_dir'
    config.set_genome_path(name, path)
    assert config.remove_genome_path(name)
    with pytest.raises(GenomeNotFoundError):
        assert config.remove_genome_path(name)


def test_list_genome_assemblies(config):
    config.set_genome_path('hg19', 'test_dir1')
    config.set_genome_path('hg18', 'test_dir2')
    assert list(config.list_genome_assemblies()) == [('hg19', 'test_dir1'),
                                                     ('hg18', 'test_dir2')]


def test_has_genome_assembly(config):
    name = 'hg19'
    path = 'test_dir'
    config.set_genome_path(name, path)
    assert config.has_genome_assembly(name)
    assert not config.has_genome_assembly('mm9')


def test_set_motif_path(config):
    name = 'motif_set'
    path = 'test_dir'
    config.set_motif_path(name, path)
    assert config._config.get('motif', name) == path


def test_get_motif_path(config):
    name = 'motif_set'
    path = 'test_dir'
    config.set_motif_path(name, path)
    assert config.get_motif_path(name) == path
    with pytest.raises(MotifSetNotFoundError):
        assert config.get_motif_path('motif_set1')


def test_remove_motif_path(config):
    name = 'motif_set'
    path = 'test_dir'
    config.set_motif_path(name, path)
    assert config.remove_motif_path(name)
    with pytest.raises(MotifSetNotFoundError):
        assert config.remove_motif_path(name)


def test_list_motif_sets(config):
    config.set_motif_path('motif_set1', 'test_dir1')
    config.set_motif_path('motif_set2', 'test_dir2')
    assert list(config.list_motif_sets()) == [('motif_set1', 'test_dir1'),
                                              ('motif_set2', 'test_dir2')]


def test_has_motif_set(config):
    name = 'motif_set'
    path = 'test_dir'
    config.set_motif_path(name, path)
    assert config.has_motif_set(name)
    assert not config.has_motif_set('motif_set1')


def test_write(config):
    config.set_genome_dir("genome_root_dir")
    config.set_motif_dir("motif_root_dir")
    config.set_genome_path('hg19', 'test_dir1')
    config.set_genome_path('mm9', 'test_dir2')
    config.set_motif_path('motif_set1', 'test_dir3')
    config.set_motif_path('motif_set2', 'test_dir4')
    config.write()
    assert os.path.isfile(config.path)
    fin = open(config.path, 'r')
    assert fin.read() == "[motifscan]\n" \
                         "genome_dir = genome_root_dir\n" \
                         "motif_dir = motif_root_dir\n\n" \
                         "[genome]\n" \
                         "hg19 = test_dir1\n" \
                         "mm9 = test_dir2\n\n" \
                         "[motif]\n" \
                         "motif_set1 = test_dir3\n" \
                         "motif_set2 = test_dir4\n\n"
