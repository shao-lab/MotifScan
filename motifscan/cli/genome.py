"""
motifscan.cli.genome
--------------------

Command line interface for the 'genome' subcommand.
"""

import logging
import os
import shutil
import sys

from motifscan.config import Config
from motifscan.exceptions import RemoteGenomeNotFoundError, \
    GenomeNotFoundError
from motifscan.genome import fasta_path_fmt, bg_freq_path_fmt, gene_path_fmt, \
    cal_bg_freq, write_bg_freq
from motifscan.genome.databases import UcscDatabase
from motifscan.io.utils import copy_file, merge_files, merge_extracted_files

logger = logging.getLogger(__name__)


def run(args, config_file=None):
    if args.list:
        config = Config(config_file)
        for name, _ in config.list_genome_assemblies():
            print(name)
        return
    if args.list_remote:
        database = UcscDatabase()
        for assembly in database.assemblies:
            print(f"{assembly.id:12}\t{database.name}\t{assembly.description}")
        return
    if args.search:
        database = UcscDatabase()
        found = False
        for assembly in database.search(args.search):
            found = True
            print(f"{assembly.id:12}\t{database.name}\t{assembly.description}")
        if not found:
            logger.info(f"No match found for {args.search!r}")
        return
    if args.install:
        install_genome(args, config_file)
        return
    if args.uninstall:
        uninstall_genome(args, config_file)
        return


def install_genome(args, config_file=None):
    config = Config(config_file)
    if config.has_genome_assembly(args.name):
        logger.error(f"Genome assembly {args.name!r} already exists!")
        sys.exit(1)

    genome_dir = os.path.abspath(
        args.output_dir or os.path.join(config.get_genome_dir(), args.name))
    logger.info(f"Installing genome assembly {args.name!r} into {genome_dir}")
    if not os.path.isdir(genome_dir):
        os.makedirs(genome_dir)
    if os.listdir(genome_dir):
        logger.error("Directory not empty! Please specify another directory "
                     "or delete files under it.")
        sys.exit(1)

    fasta_path = fasta_path_fmt.format(genome_dir, args.name)
    bg_freq_path = bg_freq_path_fmt.format(genome_dir, args.name)
    gene_path = gene_path_fmt.format(genome_dir, args.name)

    if args.remote:
        download_dir = os.path.join(genome_dir, 'downloads')
        try:
            db = UcscDatabase()
            dst_fasta = db.download_sequence(args.remote, download_dir)
            logger.debug(f"Extracting the sequence file to {fasta_path}")
            merge_extracted_files(dst_fasta, fasta_path)
            dst_gene = db.download_gene(args.remote, download_dir)
            logger.debug(f"Extracting the gene annotation file to {gene_path}")
            merge_extracted_files(dst_gene, gene_path)
            if args.clean:
                logger.debug(f"Removing the download directory {download_dir}")
                shutil.rmtree(download_dir)
        except RemoteGenomeNotFoundError as e:
            logger.error(e)
            sys.exit(1)
    else:
        logger.info("Copying the sequence file(s)")
        merge_files(args.fasta_files, fasta_path)
        logger.info("Copying the gene annotation file")
        copy_file(args.gene_file, gene_path)

    logger.info("Calculating nucleotide frequencies of the genome background")
    bg_freq = cal_bg_freq(fasta_path)
    logger.info("Writing nucleotide frequencies")
    write_bg_freq(bg_freq_path, bg_freq)

    logger.info("Updating the config file")
    config.set_genome_path(args.name, genome_dir)
    config.write()
    logger.info("Successfully installed!")


def uninstall_genome(args, config_file=None):
    try:
        config = Config(config_file)
        path = config.get_genome_path(args.uninstall)
    except GenomeNotFoundError as e:
        logger.error(e)
        sys.exit(1)

    logger.info(f"Uninstalling genome assembly {args.uninstall!r}")
    if os.path.isdir(path):
        logger.info(f"Removing files under {path}")
        try:
            shutil.rmtree(path)
        except Exception as e:
            logger.error(f"Failed to remove the genome directory: {e}")
            sys.exit(1)

    logger.info("Updating the config file")
    config.remove_genome_path(args.uninstall)
    config.write()
    logger.info("Successfully uninstalled!")
