"""
motifscan.cli.motif
-------------------

Command line interface for the 'motif' subcommand.
"""

import logging
import os
import shutil
import sys
from collections import defaultdict

import numpy as np

from motifscan.config import Config
from motifscan.exceptions import RemoteMotifPFMsNotFoundError, \
    MotifSetNotFoundError
from motifscan.genome import Genome
from motifscan.io.utils import merge_files
from motifscan.motif import pfms_path_fmt, MotifPwms, load_installed_pfms, \
    get_score_cutoffs
from motifscan.motif.databases import JasparDatabase
from motifscan.motif.cscore import c_score

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

    cutoffs_all = []
    for i in range(args.n_repeat):
        if args.n_repeat > 1:
            logger.info(
                f"Building motif score cutoffs: {i + 1} of {args.n_repeat}")
        if args.seed is not None:
            seed = args.seed + i
        else:
            seed = None
        logger.info("Random sampling background sequences")
        seqs = list(genome.random_sequences(args.n_random, max_length,
                                            args.max_n, seed))

        logger.info("Calculating background motif scores")
        matrices = [pwm.matrix.tolist() for pwm in pwms]
        sampling_scores = c_score(matrices, seqs, 3, args.n_threads)

        logger.info("Getting motif score cutoffs")
        cutoffs_all.append(get_score_cutoffs(sampling_scores))

    if args.n_repeat > 1:
        logger.info("Saving averaged motif score cutoffs")
    else:
        logger.info("Saving motif score cutoffs")

    for i, pwm in enumerate(pwms):
        cutoffs = defaultdict(list)
        for pwms_cutoffs in cutoffs_all:
            pwm_cutoffs = pwms_cutoffs[i]
            for p_value, cutoff in pwm_cutoffs.items():
                cutoffs[p_value].append(cutoff)
        for p_value in cutoffs:
            cutoff = np.mean(cutoffs[p_value])
            cutoff = np.around(cutoff, 8)
            pwm.set_cutoff(p_value=p_value, cutoff=cutoff)
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
