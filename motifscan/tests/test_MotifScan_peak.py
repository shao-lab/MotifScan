#!/bin/python
import unittest
import pandas as pd 

import sys
sys.path.insert(0,'/mnt/MAmotif/motifscan_pkg/MotifScan/lib')
import peak
import motifscan_general


class testMotifScanPeak(unittest.TestCase):
    #def setUp(self):
       #self.peak_table = pd.read_pickle('test_peak_table.pkl')
       #peaks_path = 'wgEncodeSydhTfbsHek293tZnf263UcdAlnRep1.bam.bed_peaks_100.xls'
       # self.peaks_path = 'H3K4me3_H1hesc_VS_K562_all_peak_MAvalues_3peaks.xls'
       # self.genome_path = '/home/jiawei/.MotifScan/genome/hg19'
       # self.peak_table = peak.load_peak(peaks_path,self.genome_path, 500, 'manorm')
       # self.gene_table = peak.load_ref_gene('refSeq.txt')
       
       # self.chromosome_size = pd.read_pickle('/home/jiawei/.MotifScan/genome/hg19/chromosome_size')

       # self.rand_table = pd.read_pickle('test_rnd_table.pkl')
       # self.motif_table = pd.read_pickle('test_motif_table.pkl')

    def test_generate_random_with_ref(self):
        peak.generate_random_with_ref2(self.gene_table,self.peak_table,self.genome_path,5)

    def test_merge_peak_gene(self):
        peak.merge_peak_gene(self.peak_table, self.gene_table)
    
    def test_distance2target_gene(self):
        peak_gene = peak.merge_peak_gene(self.peak_table, self.gene_table)
        print peak.distance2target_gene(peak_gene)
    
    def test_load_peak(self):
        #peaks_path = 'H3K4me3_H1hesc_VS_K562_all_peak_MAvalues_3peaks.xls'
        peaks_path = 'Ese14_H3K27ac_LICR_Rep2_peaks.xls'
        genome_path = '/home/jiawei/.MotifScan/genome/mm9'
        gene_path = '/home/jiawei/.MotifScan/gene/mm9/refSeq.txt'
        print "loading genes"
        gene_table = peak.load_ref_gene(gene_path)
        print "loading peaks"
        peak_table =  peak.load_peak(peaks_path,genome_path,0,'macs')
        print peak_table[:10]
        # print "generate randoms"
        # peak.generate_random_with_ref2(gene_table,peak_table,genome_path,1)

    def test_target_gene(self):
        peaks_path = 'Ese14_H3K27ac_LICR_Rep2_peaks_200.xls'
        genome_path = '/home/jiawei/.MotifScan/genome/mm9'
        gene_path = '/home/jiawei/.MotifScan/gene/mm9/refSeq.txt'
        print "loading genes"
        gene_table = peak.load_ref_gene(gene_path)
        print "loading peaks"
        peak_table =  peak.load_peak(peaks_path,genome_path,1000,'macs')
        peak.extract_target_gene(peak_table, gene_table)
        motifscan_general.motif_sc

def suite():
    suite = unittest.TestSuite()
    # suite.addTest(testMotifScanPeak("test_distance2target_gene"))
    suite.addTest(testMotifScanPeak("test_load_peak"))
    # suite.addTest(testMotifScanPeak("test_target_gene"))
    return suite

if __name__ == '__main__':
    unittest.main(defaultTest='suite')