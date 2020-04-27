"""
motifscan.cli.motif
-------------------

Command line interface for the 'motif' subcommand.
"""

import logging
import os
import shutil
import sys

from motifscan.config import Config
from motifscan.exceptions import RemoteMotifPFMsNotFoundError, \
    MotifSetNotFoundError
from motifscan.genome import Genome
from motifscan.io.utils import merge_files
from motifscan.motif import pfms_path_fmt, MotifPwms, load_installed_pfms
from motifscan.motif.databases import JasparDatabase
from motifscan.motif.score import motif_score

logger = logging.getLogger(__name__)


def run(args, config_file=None):
    if args.list:
        config = Config(config_file)
        for name, _ in config.list_motif_sets():
            print(name)
        return
    if args.list_remote:
        database = JasparDatabase()
        if args.database == 'jaspar_core':
            remote_sets = database.pfms_core
            for name in remote_sets:
                print(f"{name:25}\t{database.name + '_CORE'}")
        else:
            remote_sets = database.pfms_other_collections
            for name in remote_sets:
                print(f"{name:20}\t{database.name + '_Collections'}")
        return
    if args.install:
        install_motif(args, config_file)
        return
    if args.build:
        build_motif(args, config_file)
        return
    if args.uninstall:
        uninstall_motif(args, config_file)
        return


def install_motif(args, config_file=None):
    config = Config(config_file)
    if config.has_motif_set(args.name):
        logger.error(f"Motif set {args.name!r} already exists!")
        sys.exit(1)

    motif_dir = os.path.abspath(
        args.output_dir or os.path.join(config.get_motif_dir(), args.name))
    logger.info(f"Installing motif set {args.name!r} into {motif_dir}")
    if not os.path.isdir(motif_dir):
        os.makedirs(motif_dir)
    if os.listdir(motif_dir):
        logger.error("Directory not empty! Please specify another directory "
                     "or delete files under it.")
        sys.exit(1)

    pfms_path = pfms_path_fmt.format(motif_dir, args.name)

    if args.remote:
        try:
            db = JasparDatabase()
            if args.database == 'jaspar_core':
                dst_pfms = db.download_core(args.remote, motif_dir)
            else:
                dst_pfms = db.download_other_collections(args.remote,
                                                         motif_dir)
            logger.debug(
                f"Renaming downloaded file to {os.path.basename(pfms_path)}")
            shutil.move(dst_pfms, pfms_path)
        except RemoteMotifPFMsNotFoundError as e:
            logger.error(e)
            sys.exit(1)
    else:
        logger.info("Copying the PFMs file(s)")
        merge_files(args.pfm_files, pfms_path)

    logger.info("Updating the config file")
    config.set_motif_path(args.name, motif_dir)
    config.write()
    logger.info("Successfully installed!")
    if args.genome:
        build_motif(args, config_file)


def build_motif(args, config_file=None):
    if args.build:
        name = args.build
    else:
        name = args.name
    logger.info(
        f"Building motif set {name!r} for genome assembly {args.genome!r}")
    genome = Genome(args.genome)
    pfms = load_installed_pfms(name)

    logger.info("Converting motif PFMs to PWMs")
    max_length = 0
    pwms = MotifPwms(name=name, genome=genome.name)
    for pfm in pfms:
        max_length = max(max_length, pfm.length)
        pwm = pfm.to_ppm().to_pwm(genome.bg_freq)
        pwms.append(pwm)

    logger.info("Random sampling sequences")
    random_sequences = list(genome.random_sequences(args.n_random, max_length,
                                                    args.seed))
    logger.info("Scanning motifs on the sampled sequences")
    sampling_scores = []
    for idx, pwm in enumerate(pwms):
        logger.info(f"Scanning {pwm.name} [{idx + 1}/{len(pwms)}]")
        matrix = pwm.matrix.tolist()
        scores = motif_score(matrix, pwm.length, pwm.max_raw_score,
                             random_sequences)
        sampling_scores.append(scores)
    logger.info("Setting motif score cutoffs")
    pwms.set_cutoffs(sampling_scores)
    pwms.save_built_pwms()
    logger.info("Successfully built!")


def uninstall_motif(args, config_file=None):
    try:
        config = Config(config_file)
        path = config.get_motif_path(args.uninstall)
    except MotifSetNotFoundError as e:
        logger.error(e)
        sys.exit(1)

    logger.info(f"Uninstalling motif set {args.uninstall!r}")
    if os.path.isdir(path):
        logger.info(f"Removing files under {path}")
        try:
            shutil.rmtree(path)
        except Exception as e:
            logger.error(f"Failed to remove the motif directory: {e}")
            sys.exit(1)

    logger.info("Updating the config file")
    config.remove_motif_path(args.uninstall)
    config.write()
    logger.info("Successfully uninstalled!")
