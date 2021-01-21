"""
motifscan.cli.main
------------------

Main command line interface of MotifScan.
"""

import argparse
import os
import sys
from textwrap import dedent

from motifscan import __version__
from motifscan.cli import config, genome, motif, scan
from motifscan.config import user_rc_path
from motifscan.logging import setup_logger
from motifscan.region import REGION_FORMATS


def _exit(status=0, message=None):
    if message:
        print(message, file=sys.stderr)
    sys.exit(status)


def _pos_int(value):
    """Check whether a passed argument is a positive integer."""
    try:
        value_int = int(value)
        if value_int <= 0:
            raise ValueError
    except (ValueError, TypeError):
        raise argparse.ArgumentTypeError(
            f"invalid positive int value: {value!r}")
    return value_int


def _non_negative_int(value):
    """Check whether a passed argument is a non-negative integer."""
    try:
        value_int = int(value)
        if value_int < 0:
            raise ValueError
    except (ValueError, TypeError):
        raise argparse.ArgumentTypeError(
            f"invalid non-negative int value: {value!r}")
    return value_int


def _add_verbose_argument(parser):
    parser.add_argument(
        "--verbose", dest="verbose", action="store_true", default=False,
        help="Enable verbose log messages.")
    return parser


def configure_parser_main():
    """Configure the arguments parsers for MotifScan."""
    description = dedent("""
    MotifScan: A motif discovery tool to detect the occurrences of known motifs
    
    Given a set of input genomic regions, MotifScan scans the sequences to 
    detect the occurrences of known motifs. It can also perform an enrichment 
    analysis to check whether these motifs are over/under-represented compared 
    to the control regions.
    
    !!! NOTE !!!
    MotifScan requires basic data files including genome sequences and motif 
    PFMs (Position Frequency Matrices) to detect the binding sites of motifs. 
    Before scanning, users should install genome assemblies and motif sets from
    a remote database or with local prepared files via `motifscan genome` and 
    `motifscan motif` subcommands.
    
    Citation:
    Sun, H., Wang, J., Gong, Z. et al. Quantitative integration of epigenomic 
    variation and transcription factor binding using MAmotif toolkit identifies
    an important role of IRF2 as transcription activator at gene promoters. 
    Cell Discov 4, 38 (2018). https://doi.org/10.1038/s41421-018-0045-y
    """)

    epilog_msg = dedent("""
    Please run `motifscan COMMAND -h` to see the subcommand options.
     
    See also: 
      Documentation: https://motifscan.readthedocs.io
      Source code: https://github.com/shao-lab/MotifScan
      Bug reports: https://github.com/shao-lab/MotifScan/issues
    """)

    parser = argparse.ArgumentParser(
        description=description, epilog=epilog_msg,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("-v", "--version", action="version",
                        version=f"MotifScan {__version__}")

    subparsers = parser.add_subparsers(title="MotifScan Subcommands",
                                       metavar="command", dest="cmd")
    configure_parser_config(subparsers)
    configure_parser_genome(subparsers)
    configure_parser_motif(subparsers)
    configure_parser_scan(subparsers)
    return parser


def configure_parser_config(subparsers):
    """Configure the arguments parsers for 'config' subcommand."""
    help_msg = "Configure data paths for MotifScan."
    desc_msg = help_msg + dedent(f"""  

    Commands listed below enable users to change the default installation
    location of genome/motif data files and check the paths of installed 
    genome assemblies or motif sets.

    The user specific config file is located at: {user_rc_path}
    """)

    epilog_msg = dedent("""
    Examples:
    ---------    
    1) Display all values set in the config file:

        motifscan config --show

    2) Change the default installation location for genome assemblies:

        motifscan config --set-default-genome <path>

    3) Change the default installation location for motif sets:

        motifscan config --set-default-motif <path>     

    4) Get the genome path of a specific genome assembly:

        motifscan config --get-genome <genome_name>

    5) Change the motif path for a specific motif set:

        motifscan config --set-motif <motif_set> <path>
    """)

    parser = subparsers.add_parser(
        "config", description=desc_msg, help=help_msg, epilog=epilog_msg,
        formatter_class=argparse.RawDescriptionHelpFormatter)

    parser_basic = parser.add_argument_group("Basic Options")
    parser_basic.add_argument(
        "--show", dest="show", action="store_true", default=False,
        help="Show all configured values.")

    parser_default = parser.add_argument_group("Default Install Location")
    parser_default.add_argument(
        "--set-default-genome", metavar="PATH", dest="set_default_genome",
        help="Set the default installation path for genome assemblies.")
    parser_default.add_argument(
        "--set-default-motif", metavar="PATH", dest="set_default_motif",
        help="Set the default installation path for motif sets.")

    parser_genome = parser.add_argument_group("Genome Path Options")
    parser_genome.add_argument(
        "--get-genome", metavar="NAME", dest="get_genome",
        help="Get the genome path of a specific genome assembly.")
    parser_genome.add_argument(
        "--set-genome", metavar=("NAME", "PATH"), dest="set_genome", nargs=2,
        help="Set the genome path for a specific genome assembly.")
    parser_genome.add_argument(
        "--rm-genome", metavar="NAME", dest="rm_genome",
        help="Remove a specific genome assembly.")

    parser_motif = parser.add_argument_group("Motif Path Options")
    parser_motif.add_argument(
        "--get-motif", metavar="NAME", dest="get_motif",
        help="Get the motif path of a specific motif set.")
    parser_motif.add_argument(
        "--set-motif", metavar=("NAME", "PATH"), dest="set_motif", nargs=2,
        help="Set the motif path for a specific motif set.")
    parser_motif.add_argument(
        "--rm-motif", metavar="NAME", dest="rm_motif",
        help="Remove a specific motif set.")

    parser = _add_verbose_argument(parser)
    parser.set_defaults(func=config.run)


def configure_parser_genome(subparsers):
    """Configure the arguments parsers for the 'genome' subcommand."""
    help_msg = "Genome assembly commands."

    desc_msg = help_msg + dedent("""
    
    This subcommand controls the genome assemblies used by MotifScan.
    MotifScan requires a sequences FASTA file and a gene annotation file 
    (if available) for each genome assembly, users can either download them 
    from a remote database or install directly with local prepared files.
    """)

    epilog_msg = dedent("""
    Examples:
    --------- 
    1) Display installed genomes:
    
        motifscan genome --list
        
    2) Display all available genomes in a remote database:
    
        motifscan genome --list-remote
    
    3) Search genomes in a remote database by keyword (e.g. 'human'):
    
        motifscan genome --search human
    
    4) Install 'hg19' genome assembly from a remote database:
    
        motifscan genome --install -n hg19 -r hg19
                       
    5) Install 'hg19' genome assembly with local prepared files:

        motifscan genome --install -n hg19 -i <hg19.fa> -a <refGene.txt>   
        
    6) Uninstall a genome assembly:
        
        motifscan genome --uninstall <genome_name>
        
    Notes:
    ------       
    The path of newly installed genome will be automatically saved. If you 
    move the directory to another location later, please reconfigure it:
    
        motifscan config --set-genome <genome_name> <new_path>
    """)

    parser = subparsers.add_parser(
        "genome", description=desc_msg, help=help_msg, epilog=epilog_msg,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    subcommands = parser.add_argument_group("Genome Subcommands")
    subcommands = subcommands.add_mutually_exclusive_group()
    subcommands.add_argument(
        "--list", dest="list", action="store_true", default=False,
        help="Display installed genome assemblies.")
    subcommands.add_argument(
        "--list-remote", dest="list_remote", action="store_true",
        default=False, help="Display available remote genome assemblies.")
    subcommands.add_argument(
        "--search", metavar="KEYWORD", dest="search",
        help="Search for genome assemblies in a remote database.")
    subcommands.add_argument(
        "--install", dest="install", action="store_true", default=False,
        help="Install a new genome assembly.")
    subcommands.add_argument(
        "--uninstall", metavar="NAME", dest="uninstall",
        help="Uninstall a genome assembly.")
    subcommands.required = True

    parser_install = parser.add_argument_group("Install Options")
    parser_install.add_argument(
        "-n", "--name", metavar="NAME", dest="name",
        help="Name of the genome assembly to be installed.")
    parser_install.add_argument(
        "-i", metavar="FASTA", dest="fasta_files", nargs="+",
        help="Local genome sequences file(s) in FASTA format.")
    parser_install.add_argument(
        "-a", metavar="ANNOTATION", dest="gene_file",
        help="Local gene annotation (refGene.txt) file.")
    parser_install.add_argument(
        "-r", "--remote", metavar="GENOME", dest="remote",
        help="Download required data files from a remote assembly.")
    parser_install.add_argument(
        "-o", "--output-dir", metavar="DIR", dest="output_dir",
        help="Write to a given directory instead of the default directory.")

    parser_remote = parser.add_argument_group("Remote Database Options")
    parser_remote.add_argument(
        "--database", dest="database", choices=["ucsc"], default="ucsc",
        help="Which remote database is used to list/install/search genome "
             "assemblies. Default: ucsc")
    parser_remote.add_argument(
        "--clean", dest="clean", action="store_true", default=False,
        help="Clean the download directory after installation.")
    parser = _add_verbose_argument(parser)
    parser.set_defaults(func=genome.run)


def _check_args_genome(args):
    """Check the arguments of the 'genome' subcommand."""
    if args.install:
        # -n/--name must be specified
        if not args.name:
            _exit(1, "motifscan genome --install: error: argument -n/--name "
                     "is required")
        # check conflict between local model and remote mode
        if args.remote and (args.fasta_files or args.gene_file):
            _exit(1, "motifscan genome --install: error: argument -r/--remote "
                     "is not allowed with argument -i or -a")
        # -i/-a must be specified in local mode
        if not args.remote:
            if not args.fasta_files:
                _exit(1, "motifscan genome --install: error: argument -i is "
                         "required")
            if not args.gene_file:
                _exit(1, "motifscan genome --install: error: argument -a is "
                         "required")
            # check if the input files are existed
            input_files = list(args.fasta_files)
            input_files.append(args.gene_file)
            for path in input_files:
                if not os.path.isfile(path):
                    _exit(1, f"motifscan genome --install: error: file not "
                             f"found: {path}")


def configure_parser_motif(subparsers):
    """Configure the arguments parsers for the 'motif' subcommand."""
    help_msg = "Motif set (PFMs/PWMs) commands."
    desc_msg = help_msg + dedent("""
    
    MotifScan only detects the binding sites of known motifs. Before scanning, 
    the motif set should be installed and built with PFMs (Position Frequency 
    Matrices). Since different assemblies have different genome contents, it 
    is necessary to build the PFMs and get proper motif score cutoffs for every 
    genome assembly you want to scan later. 
    """)

    epilog_msg = dedent("""
    Examples:
    ---------        
    1) Display installed motif sets:
        
        motifscan motif --list
    
    2) Display all available motif sets in a remote database:
    
        motifscan motif --list-remote
    
    3) Install a motif set from a remote database and build for genome 'hg19':
    
        motifscan motif --install -n <motif_set> -r <remote_PFMs> -g hg19
         
    4） Install a motif set with local PFMs file(s) and build for genome 'mm9':

        motifscan motif --install -n <motif_set> -i <pfms.jaspar> -g mm9
    
    5) Build an installed motif set (PFMs) for additional assembly 'hg38':
    
        motifscan motif --build <motif_set> -g hg38
        
    6) Uninstall a motif set:
        
        motifscan motif --uninstall <motif_set>
        
    Notes:
    ------
    1) When installing a motif set by `--install`, you can append a `-g` option 
    to build the PFMs for the specified assembly after installation.
    
    2) The genome assembly specified by `-g` should be pre-installed by command 
    `motifscan genome --install`.
    
    3) The path of newly installed motif set will be automatically saved and 
    all the built PWMs files are stored under the directory. If you move it 
    to a new path, please reconfigure it:
    
        motifscan config --set-motif <motif_set> <new_path>
    """)

    parser = subparsers.add_parser(
        "motif", description=desc_msg, help=help_msg, epilog=epilog_msg,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    subcommands = parser.add_argument_group("Motif Subcommands")
    subcommands = subcommands.add_mutually_exclusive_group()
    subcommands.add_argument(
        "--list", dest="list", action="store_true", default=False,
        help="Display installed motif sets.")
    subcommands.add_argument(
        "--list-remote", dest="list_remote", action="store_true",
        default=False, help="Display available remote motif sets.")
    subcommands.add_argument(
        "--install", dest="install", action="store_true", default=False,
        help="Install a new motif set with PFMs.")
    subcommands.add_argument(
        "--build", metavar="NAME", dest="build", default=None,
        help="Build an installed motif set for additional genome assembly.")
    subcommands.add_argument(
        "--uninstall", metavar="NAME", dest="uninstall",
        help="Uninstall a motif set.")
    subcommands.required = True

    parser_install = parser.add_argument_group("Install Options")
    parser_install.add_argument(
        "-n", "--name", metavar="NAME", dest="name",
        help="Name of the motif set (PFMs) to be installed.")
    parser_install.add_argument(
        "-i", metavar="FILE", dest="pfm_files", nargs="+",
        help="Local motif PFMs file(s) to be installed.")
    parser_install.add_argument(
        "-r", "--remote", metavar="PFMs", dest="remote",
        help="Download a remote motif PFMs set.")
    parser_install.add_argument(
        "-o", "--output-dir", metavar="DIR", dest="output_dir",
        help="Write to a given directory instead of the default directory.")

    parser_remote = parser.add_argument_group("Remote Database Options")
    parser_remote.add_argument(
        "--database", dest="database",
        choices=["jaspar_core", "jaspar_collections"], default="jaspar_core",
        help="Which remote database is used to list/install motif set (PFMs). "
             "Default: jaspar_core")

    parser_build = parser.add_argument_group("Build Options")
    parser_build.add_argument(
        "-g", "--genome", metavar="GENOME", dest="genome",
        help="Genome assembly to build the motif set (PFMs) for.")
    parser_build.add_argument(
        "--n-random", metavar="N", dest="n_random", type=int, default=1000000,
        help="Generate N random background sequences to calculate motif score "
             "cutoffs. Default: 1,000,000")
    parser_build.add_argument(
        "--n-repeat", metavar="N", dest="n_repeat", type=_pos_int, default=1,
        help="Repeat N rounds of random sampling and use the averaged cutoff "
             "as final cutoff. Default: 1")
    parser_build.add_argument(
        "--max-n", metavar="N", dest="max_n", type=int, default=0,
        help="The maximal number of `N` base allowed in each random sampled "
             "sequence. Default: 0")
    parser_build.add_argument(
        "--seed", metavar="SEED", dest="seed", type=int, default=None,
        help="Random seed used to generate background sequences.")

    parser_threads = parser.add_argument_group("Threads Options")
    parser_threads.add_argument(
        "-t", "--threads", metavar="N", dest="n_threads", type=int, default=1,
        help="Number of processes used to run in parallel.")

    parser = _add_verbose_argument(parser)
    parser.set_defaults(func=motif.run)


def _check_args_motif(args):
    """Check the arguments of the 'motif' subcommand."""
    if args.install:
        # -n/--name must be specified
        if not args.name:
            _exit(1, "motifscan motif --install: error: argument -n/--name "
                     "is required")
        # check conflict between local model and remote mode
        if args.remote and args.pfm_files:
            _exit(1, "motifscan motif --install: error: argument -r/--remote "
                     "is not allowed with argument -i")
        # -i must be specified in local mode
        if not args.remote:
            if not args.pfm_files:
                _exit(1, "motifscan motif --install: error: argument -i is "
                         "required")
            # check the input files are existed
            for path in args.pfm_files:
                if not os.path.isfile(path):
                    _exit(1, f"motifscan motif --install: error: file not "
                             f"found: {path}")
    if args.build:
        # -g/--genome must be specified
        if not args.genome:
            _exit(1, "motifscan motif --build: error: argument -g/--genome "
                     "is required")


def configure_parser_scan(subparsers):
    """Configure the arguments parsers for the 'scan' subcommand."""
    help_msg = "Scan input regions to detect motif occurrences."
    desc_msg = help_msg + dedent("""
    
    This main command invokes to scan the sequences of user specified input 
    genomic regions and detect the occurrences for a set of known motifs. 
    After scanning the input regions, an optional motif enrichment analysis 
    is performed to check whether these motifs are over/under-represented 
    compared to control regions (can be random generated or user specified).
    """)

    epilog_msg = dedent("""
    Examples:
    ---------
    1) Scan input regions for a set of known motifs under 'hg19' genome:

        motifscan scan -i regions.bed -m <motif_set> -g hg19 -o <path>
        
    2) Test motif enrichment compared to user-specified control regions:

        motifscan scan -i regions.bed -c control.bed -m <motif_set> -g hg19 -o <path>
    
    3) Only scan input regions located at promoters:

        motifscan scan -i regions.bed -m <motif_set> -g hg19 --loc promoter -o <path>
        
    4) Scan whole input regions rather than fixed-size windows:

        motifscan scan -i regions.bed -m <motif_set> -g hg19 -w 0 -o <path>

    5) Report the positions and distributions of detected motif sites:

        motifscan scan -i regions.bed -m <motif_set> -g hg19 --site --plot -o <path>

    """)

    parser = subparsers.add_parser(
        "scan", description=desc_msg, help=help_msg, epilog=epilog_msg,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser_input = parser.add_argument_group("Input Options")
    parser_input.add_argument(
        "-i", metavar="FILE", dest="input_file", required=True,
        help="Input genomic regions (peaks) to be scanned.")
    parser_input.add_argument(
        "-f", dest="input_format", choices=REGION_FORMATS, default="bed",
        help="Format of the input file. Default: bed")
    parser_input.add_argument(
        "-m", "--motif", metavar="NAME", dest="motif", required=True,
        help="Motif set name to scan for.")
    parser_input.add_argument(
        "-g", "--genome", metavar="GENOME", dest="genome", required=True,
        help="Genome assembly name.")

    parser_advance = parser.add_argument_group("Scanning Options")
    parser_advance.add_argument(
        "-p", dest="p_value", default="1e-4",
        choices=["1e-2", "1e-3", "1e-4", "1e-5", "1e-6"],
        help="P value cutoff for motif scores. Default: 1e-4")
    parser_advance.add_argument(
        "--loc", dest="location", choices=["promoter", "distal"], default=None,
        help="If specified, only scan promoter or distal regions.")
    parser_advance.add_argument(
        "--upstream", metavar="DISTANCE", dest="upstream",
        type=_pos_int, default=4000,
        help="TSS upstream distance for promoters. Default: 4000")
    parser_advance.add_argument(
        "--downstream", metavar="DISTANCE", dest="downstream",
        type=_pos_int, default=2000,
        help="TSS downstream distance for promoters. Default: 2000")
    parser_advance.add_argument(
        "-w", "--window-size", metavar="LENGTH", dest="window_size",
        type=_non_negative_int, default=1000,
        help="Window size for scanning. In most cases, motifs occur closely "
             "around the centers or summits of genomic peaks. Scanning a "
             "fixed-size window is often sufficient to detect motif sites and "
             "unbiased for the enrichment analysis. If set to 0, the whole "
             "input regions are included for scanning. Default: 1000")
    parser_advance.add_argument(
        "--strand", dest="strand", choices=['both', '+', '-'], default="both",
        help="Enable strand-specific scanning, defaults to scan both strands.")

    parser_enrich = parser.add_argument_group("Enrichment Analysis Options")
    parser_enrich.add_argument(
        "--no-enrich", dest="no_enrich", action="store_true", default=False,
        help="Disable the enrichment analysis.")
    parser_enrich.add_argument(
        "--n-random", metavar="N", dest="n_random",
        type=_non_negative_int, default=5,
        help="Generate N random control regions for each input region. "
             "Default: 5")
    parser_enrich.add_argument(
        "--seed", metavar="SEED", dest="seed", type=int, default=None,
        help="Random seed used to generate control regions.")
    parser_enrich.add_argument(
        "-c", metavar="FILE", dest="control_file",
        help="Use custom control regions for the enrichment analysis.")
    parser_enrich.add_argument(
        "--cf", dest="control_format", choices=REGION_FORMATS, default="bed",
        help="Format of the control file. Default: bed")

    parser_threads = parser.add_argument_group("Threads Options")
    parser_threads.add_argument(
        "-t", "--threads", metavar="N", dest="n_threads", type=int, default=1,
        help="Number of processes used to run in parallel.")

    parser_output = parser.add_argument_group("Output Options")
    parser_output.add_argument(
        "-o", "--output-dir", metavar="DIR", dest="output_dir", required=True,
        help="Directory to write output files.")
    parser_output.add_argument(
        "--site", dest="report_site", action="store_true", default=False,
        help="If set, report the position for each detected motif site.")
    parser_output.add_argument(
        "--plot", dest="plot_dist", action="store_true", default=False,
        help="If set, plot the distributions of detected motif sites.")

    parser = _add_verbose_argument(parser)
    parser.set_defaults(func=scan.run)


def main():
    parser = configure_parser_main()
    args = parser.parse_args()
    if args.cmd == 'genome':
        _check_args_genome(args)
    elif args.cmd == 'motif':
        _check_args_motif(args)
    setup_logger(args.verbose)
    args.func(args)


if __name__ == '__main__':
    main()
