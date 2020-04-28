"""
motifscan.cli.scan
------------------

Command line interface for the 'scan' subcommand.
"""

import logging
import sys

from motifscan import __version__
from motifscan.genome import Genome
from motifscan.io import write_sites_table, write_sites_bed, write_enrich_table
from motifscan.motif import load_built_pwms
from motifscan.plot import plot_motif_sites_dist, plot_motif_sites_enrich
from motifscan.region import load_motifscan_regions
from motifscan.region.utils import subset_by_location, generate_control_regions
from motifscan.scanner import Scanner
from motifscan.stats import motif_enrichment

logger = logging.getLogger(__name__)


def run(args):
    logger.info(f"Running MotifScan {__version__}")
    logger.info("===== Loading data =====")
    genome = Genome(name=args.genome)
    pwms = load_built_pwms(name=args.motif, genome=args.genome)
    regions = load_motifscan_regions(path=args.input_file,
                                     format=args.input_format)
    if args.location is not None:
        logger.info(f"Extracting input regions located at {args.location}")
        if genome.genes is None:
            logger.error(f"Unable to extract without gene annotations.")
            sys.exit(1)
        regions = subset_by_location(
            regions=regions, genes=genome.genes, location=args.location,
            upstream=args.upstream, downstream=args.downstream)
        logger.info(f"Extracted {len(regions)} {args.location} regions")

    logger.info("===== Scanning motifs =====")
    logger.info("Fetching the sequences of input regions")
    scanner = Scanner(
        genome=genome, regions=regions, window_size=args.window_size,
        strand=args.strand, p_value=args.p_value, remove_dup=True,
        n_threads=args.n_threads)
    logger.info("Scanning motifs...")
    motif_sites = scanner.scan_motifs(pwms=pwms)

    logger.info("Saving the result tables")
    write_sites_table(output_dir=args.output_dir, pwms=pwms, regions=regions,
                      motif_sites=motif_sites)
    if args.report_site:
        logger.info("Saving the coordinates of detected motif sites")
        write_sites_bed(output_dir=args.output_dir, pwms=pwms, regions=regions,
                        motif_sites=motif_sites)

    if not args.no_enrich:
        logger.info("===== Motif Enrichment =====")
        if args.control_file:
            logger.info("Loading user specified control regions")
            control_regions = load_motifscan_regions(
                path=args.control_file, format=args.control_format)
            if args.location is not None:
                logger.info(
                    f"Extracting control regions located at {args.location}")
                control_regions = subset_by_location(
                    regions=control_regions, genes=genome.genes,
                    location=args.location,
                    upstream=args.upstream, downstream=args.downstream)
                logger.info(
                    f"Extracted {len(control_regions)} {args.location} "
                    f"control regions")
        else:
            logger.info(f"Generating random control regions")
            control_regions = generate_control_regions(
                n_random=args.n_random, regions=regions,
                chrom_size=genome.chrom_sizes, genes=genome.genes,
                random_seed=args.seed)
        logger.info("Fetching the sequences of control regions")
        scanner_control = Scanner(
            genome=genome, regions=control_regions,
            window_size=args.window_size, strand=args.strand,
            p_value=args.p_value, remove_dup=True, n_threads=args.n_threads)
        logger.info("Scanning motifs...")
        motif_sites_control = scanner_control.scan_motifs(pwms=pwms)

        logger.info("Performing motif enrichment analysis")
        enrichment_results = motif_enrichment(
            pwms=pwms, motif_sites=motif_sites,
            motif_sites_control=motif_sites_control)
        logger.info("Saving the motif enrichment table")
        write_enrich_table(output_dir=args.output_dir,
                           enrichment_results=enrichment_results)

    if args.plot_dist:
        logger.info("Plotting the distributions of detected motif sites")
        plot_motif_sites_dist(
            output_dir=args.output_dir, regions=regions, pwms=pwms,
            motif_sites=motif_sites, window_size=args.window_size)
        if not args.no_enrich:
            logger.info("Plotting the enrichment of detected motif sites")
            plot_motif_sites_enrich(
                output_dir=args.output_dir, regions=regions, pwms=pwms,
                motif_sites=motif_sites,
                motif_sites_control=motif_sites_control)

    logger.info("===== MotifScan Finished =====")
