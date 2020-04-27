"""
motifscan.genome.annotation
---------------------------

Module for gene annotations.
"""

import logging
from collections import defaultdict

logger = logging.getLogger(__name__)


class Gene:
    """Class for individual gene (transcript)."""

    def __init__(self, chrom, tss, strand, name=None):
        self.chrom = chrom
        self.tss = int(tss)
        if strand not in ['+', '-']:
            raise ValueError(f"invalid strand option: {strand!r}")
        self.strand = strand
        self.name = name

    def promoter(self, upstream=2000, downstream=2000):
        if self.strand == '+':
            return [self.tss - upstream, self.tss + downstream]
        else:
            return [self.tss - downstream, self.tss + upstream]


class Genes:
    """Class for a set of genes."""

    def __init__(self, path):
        self.path = path
        self._genes = defaultdict(list)
        self.read_genes()

    def __len__(self):
        return sum(len(self._genes[chrom]) for chrom in self._genes)

    def fetch(self, chrom):
        if chrom in self._genes:
            return self._genes[chrom]
        else:
            return []

    def read_genes(self):
        logger.debug(f"Loading genes from {self.path}")
        parser = RefGeneTxtParser(self.path)
        for gene in parser.parse():
            self._genes[gene.chrom].append(gene)
        logger.debug(f"Loaded {len(self)} genes")


class RefGeneTxtParser:
    """Parser for refGene in txt format."""

    def __init__(self, path):
        self.path = path

    def parse(self):
        with open(self.path, 'r') as fin:
            for line in fin:
                line = line.strip()
                fields = line.split()
                name = fields[1]
                chrom = fields[2]
                strand = fields[3]
                if strand == '+':
                    tss = int(fields[4])
                elif strand == '-':
                    tss = int(fields[5])
                else:
                    raise ValueError(
                        f"Invalid strand {strand!r} detected at line: {line}")
                yield Gene(chrom=chrom, tss=tss, strand=strand, name=name)


def read_gene_annotation(path):
    return Genes(path)
