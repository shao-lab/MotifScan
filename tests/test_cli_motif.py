import os

import pytest

from motifscan.cli.main import configure_parser_main
from motifscan.cli.motif import run
from motifscan.config import Config


parser = configure_parser_main()


def test_cli_motif_list(tmp_dir, capsys):
    config_file = os.path.join(tmp_dir, "test_cli_motif.motifscanrc")
    config = Config(config_file)
    config.set_motif_path("motif_set1", "/path/to/motif1")
    config.set_motif_path("motif_set2", "/path/to/motif2")
    config.write()

    args = parser.parse_args(["motif", "--list"])
    run(args=args, config_file=config_file)
    captured = capsys.readouterr()
    assert captured.out == "motif_set1\nmotif_set2\n"


def test_cli_motif_list_remote(capsys):
    args = parser.parse_args(["motif", "--list-remote"])
    run(args=args)
    captured = capsys.readouterr()
    assert captured.out
    args = parser.parse_args(
        ["motif", "--list-remote", "--database", "jaspar_collections"])
    run(args=args)
    captured = capsys.readouterr()
    assert captured.out


def test_cli_motif_install(motif_root, tmp_dir):
    config_file = os.path.join(tmp_dir, "test_cli_motif.motifscanrc")
    config = Config(config_file)
    config.set_motif_dir(tmp_dir)
    config.write()

    pfms_path = os.path.join(motif_root, "test", "test_pfms.jaspar")
    args = parser.parse_args(
        ["motif", "--install", "-n", "test_motif", "-i", pfms_path])
    run(args=args, config_file=config_file)

    motif_path = os.path.join(tmp_dir, "test_motif")
    assert os.path.isfile(os.path.join(motif_path, "test_motif_pfms.jaspar"))

    config = Config(config_file)
    assert config.has_motif_set("test_motif")
    assert config.get_motif_path("test_motif") == motif_path


def test_cli_motif_uninstall(motif_root, tmp_dir):
    config_file = os.path.join(tmp_dir, "test_cli_motif.motifscanrc")
    config = Config(config_file)
    config.set_motif_dir(tmp_dir)
    config.write()

    args = parser.parse_args(["motif", "--uninstall", "test_motif"])
    run(args=args, config_file=config_file)

    config = Config(config_file)
    assert not config.has_motif_set("test_motif")

    motif_path = os.path.join(tmp_dir, "test_motif")
    assert not os.path.isdir(motif_path)

    args = parser.parse_args(["motif", "--uninstall", "test_motif1"])
    with pytest.raises(SystemExit):
        run(args=args, config_file=config_file)
