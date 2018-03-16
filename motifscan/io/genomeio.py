"""Module script for genome IO."""

import numpy as np
import pandas as pd


def dump_genome(bg_freq, chr_size, output_dir):
    """Dump the pre-compiled genome file to pickle.

        Args:
            bg_freq (pd.Series): Genome background frequency of nucleotides.
            chr_size (pd.Series): Chromosome size.
            output_dir (str): Path of pre-compiled genome directory.

        """
    bg_freq.to_pickle("{}/background".format(output_dir))
    bg_freq.to_csv("{}/background.txt".format(output_dir), sep='\t', header=False, index=True)
    chr_size.to_pickle("{}/chromosome_size".format(output_dir))
    chr_size.to_csv("{}/chromosome_size.txt".format(output_dir), sep='\t', header=False, index=True)
    return


def load_genome(genome_dir):
    """Load pre-compiled genome file from pickle.
    
    Args:
        genome_dir (str): Path of pre-compiled genome directory. 
    
    Returns:
        bg_freq (pd.Series): Genome background frequency of nucleotides.
        chr_size (pd.Series): Chromosome size.
    
    """
    bg_freq = pd.read_pickle("{}/background".format(genome_dir))
    chr_size = pd.read_pickle("{}/chromosome_size".format(genome_dir))
    return bg_freq, chr_size


def split_genome(genome_file, output_dir):
    """Split the genome sequences by chromosome.

    Args: 
        genome_file (str): Path of genome sequences input file.
        output_dir (str): Path of pre-compiled genome directory to write the sequences under.

    """
    chrs = []
    with open(genome_file, 'r') as fin:
        tmp_chr = None
        seq_buffer = ''
        for line in fin:
            if line[0] == '>':  # chromosome id line
                if tmp_chr is not None:
                    write_sequence(output_dir, tmp_chr, seq_buffer)
                tmp_chr = line.strip()[1:]
                chrs.append(tmp_chr)
                seq_buffer = line
            else:  # sequence line
                seq_buffer += line
        if tmp_chr is not None:
            write_sequence(output_dir, tmp_chr, seq_buffer)
    return


def write_sequence(output_dir, chr_name, sequence):
    """Write genome sequence of single chromosome with specific format.

    Args:
        output_dir (str): Path of pre-compiled genome directory to write the sequences under.
        chr_name (str): Chromosome name.
        sequence (str): Sequences of the chromosome.

    """
    with open("{0}/{1}".format(output_dir, chr_name), 'w') as fout:
        header_idx = sequence.find('\n') + 1
        header_line = sequence[:header_idx]
        fout.write(header_line)
        tmp_seq = sequence[header_idx:].replace('\n', '')
        tmp_len = len(tmp_seq)
        line_num = int(np.ceil(tmp_len / 50.0))
        for i in xrange(line_num):
            start_idx = 50 * i
            end_idx = 50 * (i + 1)
            fout.write("{}\n".format(tmp_seq[start_idx:end_idx]))
        fout.close()
    return
