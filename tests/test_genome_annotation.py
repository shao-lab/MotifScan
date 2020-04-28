import os

import pytest

from motifscan.genome.annotation import Gene, Genes, RefGeneTxtParser, \
    read_gene_annotation


def test_gene_init():
    gene = Gene(chrom='chr1', tss=10000, strand='+', name='gene1')
    assert gene.chrom == 'chr1'
    assert gene.tss == 10000
    assert gene.strand == '+'
    assert gene.name == 'gene1'
    assert gene.promoter(upstream=4000, downstream=2000) == [6000, 12000]
    gene = Gene(chrom='chr1', tss=10000, strand='-', name='gene2')
    assert gene.promoter(upstream=4000, downstream=2000) == [8000, 14000]
    with pytest.raises(ValueError):
        Gene(chrom='chr1', tss=10000, strand='.', name='bad')


def test_genes_init(genome_root):
    gene_path = os.path.join(genome_root, 'test', 'test_gene_annotation.txt')
    genes = Genes(gene_path)
    assert genes.path == gene_path
    assert len(genes) == 5
    genes_chr1 = genes.fetch('chr1')
    assert len(genes_chr1) == 4
    assert genes.fetch('chrX') == []


def test_ref_gene_txt_parser(genome_root):
    gene_path = os.path.join(genome_root, 'test', 'test_gene_annotation.txt')
    parser = RefGeneTxtParser(gene_path)
    assert parser.path == gene_path
    genes = list(parser.parse())
    assert len(genes) == 5
    assert genes[0].chrom == 'chr1'
    assert genes[0].tss == 11868
    assert genes[0].strand == '+'
    assert genes[0].name == 'NR_148357'
    assert genes[1].chrom == 'chr1'
    assert genes[1].tss == 11873
    assert genes[1].strand == '+'
    assert genes[1].name == 'NR_046018'
    assert genes[2].chrom == 'chr22'
    assert genes[2].tss == 24666798
    assert genes[2].strand == '+'
    assert genes[2].name == 'NM_015330'
    assert genes[3].chrom == 'chr1'
    assert genes[3].tss == 17436
    assert genes[3].strand == '-'
    assert genes[3].name == 'NR_106918'
    assert genes[4].chrom == 'chr1'
    assert genes[4].tss == 17436
    assert genes[4].strand == '-'
    assert genes[4].name == 'NR_107062'
    with pytest.raises(ValueError):
        gene_path = os.path.join(genome_root, 'bad',
                                 'test_gene_annotation_bad.txt')
        parser = RefGeneTxtParser(gene_path)
        list(parser.parse())


def test_read_gene_annotation(genome_root):
    gene_path = os.path.join(genome_root, 'test', 'test_gene_annotation.txt')
    genes = read_gene_annotation(gene_path)
    assert genes.path == gene_path
    assert len(genes) == 5
    genes_chr1 = genes.fetch('chr1')
    assert len(genes_chr1) == 4
    assert genes.fetch('chrX') == []
