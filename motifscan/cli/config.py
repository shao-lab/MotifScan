"""
motifscan.cli.config
--------------------

Command line interface for the 'config' subcommand.
"""

import logging
import os
import sys

from motifscan.config import Config
from motifscan.exceptions import GenomeNotFoundError, MotifSetNotFoundError

logger = logging.getLogger(__name__)


def run(args, config_file=None):
    config = Config(config_file)
    modified = False

    if args.show:
        print("[motifscan]")
        print(f"genome_dir: {config.get_genome_dir()}")
        print(f"motif_dir: {config.get_motif_dir()}")
        print("\n[genome]")
        for name, path in config.list_genome_assemblies():
            print(f"{name}: {path}")
        print("\n[motif]")
        for name, path in config.list_motif_sets():
            print(f"{name}: {path}")
        return

    if args.set_default_genome:
        logger.debug(
            "Setting the default installation path for genome assemblies")
        config.set_genome_dir(os.path.abspath(args.set_default_genome))
        modified = True
    if args.set_default_motif:
        logger.debug("Setting the default installation path for motif sets")
        config.set_motif_dir(os.path.abspath(args.set_default_motif))
        modified = True

    if args.get_genome:
        logger.debug(f"Getting the genome path of {args.get_genome!r}")
        try:
            print(config.get_genome_path(args.get_genome))
        except GenomeNotFoundError as e:
            logger.error(e)
            sys.exit(1)
    if args.set_genome:
        name = args.set_genome[0]
        path = os.path.abspath(args.set_genome[1])
        logger.debug(f"Setting the genome path for {name!r}: {path}")
        config.set_genome_path(name, path)
        modified = True
    if args.rm_genome:
        logger.debug(f"Removing the genome path for {args.rm_genome!r}")
        try:
            config.remove_genome_path(args.rm_genome)
            modified = True
        except GenomeNotFoundError as e:
            logger.error(e)
            sys.exit(1)

    if args.get_motif:
        logger.debug(f"Getting the motif path of {args.get_motif!r}")
        try:
            print(config.get_motif_path(args.get_motif))
        except MotifSetNotFoundError as e:
            logger.error(e)
            sys.exit(1)
    if args.set_motif:
        name = args.set_motif[0]
        path = os.path.abspath(args.set_motif[1])
        logger.debug(f"Setting the motif path for {name!r}: {path}")
        config.set_motif_path(name, path)
        modified = True
    if args.rm_motif:
        logger.debug(f"Removing the motif path for {args.rm_motif!r}")
        try:
            config.remove_motif_path(args.rm_motif)
            modified = True
        except MotifSetNotFoundError as e:
            logger.error(e)
            sys.exit(1)

    if modified:
        logger.debug(f"Updating the config file: {config.path}")
        config.write()
        logger.debug("Done")
