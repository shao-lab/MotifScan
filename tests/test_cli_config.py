import os

import pytest

from motifscan.cli.config import run
from motifscan.cli.main import configure_parser_main
from motifscan.config import Config
from motifscan.exceptions import GenomeNotFoundError, MotifSetNotFoundError

parser = configure_parser_main()


def test_cli_config_show(tmp_dir, capsys):
    config_file = os.path.join(tmp_dir, "test_cli_config.motifscanrc")
    config = Config(config_file)
    config.set_genome_dir("/path/to/genome/root")
    config.set_motif_dir("/path/to/motif/root")
    config.set_genome_path("hg19", "/path/to/genome")
    config.set_motif_path("motif_set", "/path/to/motif")
    config.write()

    args = parser.parse_args(["config", "--show"])
    run(args=args, config_file=config_file)
    captured = capsys.readouterr()
    assert captured.out == "[motifscan]\n" \
                           "genome_dir: /path/to/genome/root\n" \
                           "motif_dir: /path/to/motif/root\n\n" \
                           "[genome]\n" \
                           "hg19: /path/to/genome\n\n" \
                           "[motif]\n" \
                           "motif_set: /path/to/motif\n"


def test_cli_config_set_default_genome(tmp_dir):
    args = parser.parse_args(
        ["config", "--set-default-genome", "/path/to/genome/root"])
    config_file = os.path.join(tmp_dir, "test_cli_config.motifscanrc")
    run(args=args, config_file=config_file)

    config = Config(config_file)
    assert config.get_genome_dir() == "/path/to/genome/root"


def test_cli_config_set_default_motif(tmp_dir):
    args = parser.parse_args(
        ["config", "--set-default-motif", "/path/to/motif/root"])
    config_file = os.path.join(tmp_dir, "test_cli_config.motifscanrc")
    run(args=args, config_file=config_file)

    config = Config(config_file)
    assert config.get_motif_dir() == "/path/to/motif/root"


def test_cli_config_set_genome(tmp_dir):
    args = parser.parse_args(
        ["config", "--set-genome", "hg19", "/path/to/genome"])
    config_file = os.path.join(tmp_dir, 'test_cli_config.motifscanrc')
    run(args=args, config_file=config_file)

    config = Config(config_file)
    assert config.has_genome_assembly("hg19")
    assert config.get_genome_path("hg19") == "/path/to/genome"


def test_cli_config_get_genome(tmp_dir, capsys):
    config_file = os.path.join(tmp_dir, "test_cli_config.motifscanrc")
    config = Config(config_file)
    config.set_genome_path("hg19", "/path/to/genome")
    config.write()

    args = parser.parse_args(["config", "--get-genome", "hg19"])
    run(args=args, config_file=config_file)
    captured = capsys.readouterr()
    assert captured.out == "/path/to/genome\n"

    with pytest.raises(SystemExit):
        args = parser.parse_args(["config", "--get-genome", "hg38"])
        run(args=args, config_file=config_file)


def test_cli_config_rm_genome(tmp_dir):
    config_file = os.path.join(tmp_dir, "test_cli_config.motifscanrc")
    config = Config(config_file)
    config.set_genome_path("hg19", "/path/to/genome")
    config.write()

    args = parser.parse_args(["config", "--rm-genome", "hg19"])
    run(args=args, config_file=config_file)
    config = Config(config_file)
    assert not config.has_genome_assembly("hg19")

    with pytest.raises(GenomeNotFoundError):
        config.get_genome_path("hg19")
    with pytest.raises(SystemExit):
        run(args=args, config_file=config_file)


def test_cli_config_set_motif(tmp_dir):
    args = parser.parse_args(
        ["config", "--set-motif", "motif_set", "/path/to/motif"])
    config_file = os.path.join(tmp_dir, 'test_cli_config.motifscanrc')
    run(args=args, config_file=config_file)

    config = Config(config_file)
    assert config.has_motif_set("motif_set")
    assert config.get_motif_path("motif_set") == "/path/to/motif"


def test_cli_config_get_motif(tmp_dir, capsys):
    config_file = os.path.join(tmp_dir, "test_cli_config.motifscanrc")
    config = Config(config_file)
    config.set_motif_path("motif_set", "/path/to/motif")
    config.write()

    args = parser.parse_args(["config", "--get-motif", "motif_set"])
    run(args=args, config_file=config_file)
    captured = capsys.readouterr()
    assert captured.out == "/path/to/motif\n"

    with pytest.raises(SystemExit):
        args = parser.parse_args(["config", "--get-motif", "motif_set1"])
        run(args=args, config_file=config_file)


def test_cli_config_rm_motif(tmp_dir):
    config_file = os.path.join(tmp_dir, "test_cli_config.motifscanrc")
    config = Config(config_file)
    config.set_motif_path("motif_set", "/path/to/motif")
    config.write()

    args = parser.parse_args(["config", "--rm-motif", "motif_set"])
    run(args=args, config_file=config_file)
    config = Config(config_file)
    assert not config.has_motif_set("motif_set")

    with pytest.raises(MotifSetNotFoundError):
        config.get_motif_path("motif_set")
    with pytest.raises(SystemExit):
        run(args=args, config_file=config_file)
