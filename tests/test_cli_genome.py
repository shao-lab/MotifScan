import os

import pytest

from motifscan.cli.genome import run
from motifscan.cli.main import configure_parser_main
from motifscan.config import Config
from motifscan.genome import Genome

parser = configure_parser_main()


def test_cli_genome_list(tmp_dir, capsys):
    config_file = os.path.join(tmp_dir, "test_cli_genome.motifscanrc")
    config = Config(config_file)
    config.set_genome_path("hg19", "/path/to/genome1")
    config.set_genome_path("hg38", "/path/to/genome2")
    config.write()

    args = parser.parse_args(["genome", "--list"])
    run(args=args, config_file=config_file)
    captured = capsys.readouterr()
    assert captured.out == "hg19\nhg38\n"


def test_cli_genome_list_remote(capsys):
    args = parser.parse_args(["genome", "--list-remote"])
    run(args=args)
    captured = capsys.readouterr()
    assert captured.out


def test_cli_genome_search(capsys):
    args = parser.parse_args(["genome", "--search", "human"])
    run(args=args)
    captured = capsys.readouterr()
    assert captured.out


def test_cli_genome_install(genome_root, tmp_dir):
    config_file = os.path.join(tmp_dir, "test_cli_genome.motifscanrc")
    config = Config(config_file)
    config.set_genome_dir(tmp_dir)
    config.write()

    fasta_path = os.path.join(genome_root, "test", "test.fa")
    gene_path = os.path.join(genome_root, "test", "test_gene_annotation.txt")
    args = parser.parse_args(
        ["genome", "--install", "-n", "test_genome", "-i", fasta_path, "-a",
         gene_path])
    run(args=args, config_file=config_file)

    genome_path = os.path.join(tmp_dir, "test_genome")
    genome = Genome(name="test_genome", path=genome_path)
    assert genome.fetch_sequence("chr1", 0, 10) == "AaTtCcGgNn"
    assert genome.genes

    config = Config(config_file)
    assert config.has_genome_assembly("test_genome")
    assert config.get_genome_path("test_genome") == genome_path


def test_cli_genome_uninstall(genome_root, tmp_dir):
    config_file = os.path.join(tmp_dir, "test_cli_genome.motifscanrc")
    config = Config(config_file)
    config.set_genome_dir(tmp_dir)
    config.write()

    args = parser.parse_args(["genome", "--uninstall", "test_genome"])
    run(args=args, config_file=config_file)

    config = Config(config_file)
    assert not config.has_genome_assembly("test_genome")

    genome_path = os.path.join(tmp_dir, "test_genome")
    assert not os.path.isdir(genome_path)

    args = parser.parse_args(["genome", "--uninstall", "test_genome1"])
    with pytest.raises(SystemExit):
        run(args=args, config_file=config_file)
