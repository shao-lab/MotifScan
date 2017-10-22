import unittest
import pandas as pd
import numpy as np  
import sys
sys.path.insert(0,'/mnt/MAmotif/motifscan_pkg/MotifScan/lib')
import core
import peak
import motifscan_general
import motif
import genome
from scipy import signal
from ctypes import * 
import copy


class testMotifScanCore(unittest.TestCase):

    def test_tarsite_fc_Plot(self):
        self.peak_table = pd.read_pickle('test_peak_table.pkl')
        self.rand_table = pd.read_pickle('test_rnd_table.pkl')
        self.motif_table = pd.read_pickle('test_motif_table.pkl')
        core.tarnum_and_tarsite_distribution(self.peak_table, self.rand_table, self.motif_table, '.', bin_size=5, region_radius=500)
    
    def test_sliding_score(self):
        s = range(10)
        w = [0.1, 0.2, 0.3]
        print core.sliding_score(w,s) 
    
    def test_lfilter_speed(self):
        s = range(1000000)
        w = [0.3, 0.2, 0.1]
        signal.lfilter(w,1,s)
  
    def test_sliding_score_mat(self):
        score_c = CDLL("../lib/score.so")
        W = np.array([1,2,3,4,5,6,1,2,3,4,5,6])
        S = np.array([1,2,3,4,1,2,3,4,1,2,3,4,1,2,3,4])
        # tmp = c_double*12
        # W = tmp(1,2,3,4,5,6,1,2,3,4,5,6)
        # tmp = c_double*16
        # S = tmp(1,2,3,4,1,2,3,4,1,2,3,4,1,2,3,4)
        W = W.ctypes.data_as(POINTER(c_double))
        S = S.ctypes.data_as(POINTER(c_double))
        score_c.sliding_score_mat.restype = c_double
        print score_c.sliding_score_mat(W,c_int(3),S,c_int(4))

    def test_c_sliding_score_mat(self):
        class MAT(Structure):
            _fields_ = [("n", c_int),
                       ("a_arr",POINTER(c_double)),
                       ("c_arr",POINTER(c_double)),
                       ("g_arr",POINTER(c_double)),
                       ("t_arr",POINTER(c_double))]
        score_c = CDLL("../lib/score.so")
        score_c.sliding_score_mat.restype = POINTER(c_double)
        score_c.sliding_score_mat.argtypes = [POINTER(MAT),POINTER(MAT)]
        a_s = (c_double*6)(2,3,4,1,2,3)
        c_s = (c_double*6)(2,3,4,5,3,1)
        g_s = (c_double*6)(3,5,6,1,2,3)
        t_s = (c_double*6)(3,5,1,2,2,3)

        a = (c_double*3)(0.2,0.2,0.2)
        c = (c_double*3)(0.1,0.3,0.2)
        g = (c_double*3)(0.4,0.2,0.5)
        t = (c_double*3)(0.3,0.3,0.1)
       
        mat = MAT(3,a,c,g,t)
        smatrix = MAT(6,a_s,c_s,g_s,t_s)
        score_ptr = pointer(score_c.sliding_score_mat(byref(mat),byref(smatrix)))
        print score_ptr.contents[0]


    def test_c_testMAT(self):
        class MAT(Structure):
            _fields_ = [("n", c_int),
                       ("a_arr",POINTER(c_double)),
                       ("c_arr",POINTER(c_double)),
                       ("g_arr",POINTER(c_double)),
                       ("t_arr",POINTER(c_double))]
        score_c = CDLL("../lib/score.so")
        score_c.testMAT.restype = POINTER(MAT)
       
        score_ptr = score_c.testMAT()
        print score_ptr.contents.a_arr

    def test_c_testMOTIF_RES(self):
        class MOTIF_RES(Structure):
            _fields_ = [("tarnum", c_int),
                       ("ratio", c_double),
                       ("tarsite",POINTER(c_int)),
                       ("tarratio", POINTER(c_double))]
        score_c = CDLL("../lib/score.so")
        score_c.testMOTIF_RES.restype = POINTER(MOTIF_RES)
       
        score_ptr = score_c.testMOTIF_RES()
        print score_ptr.contents.tarsite[0]
        print score_ptr.contents.tarsite[1]
        print score_ptr.contents.tarsite[2]
        for i in np.arange(10):
            print score_ptr.contents.tarratio[i]


    def test_motif_scan_core(self):
        max_score = 0.8
        score_cutoff = -1
        B = [0.28,0.22,0.22,0.28]
        # for i in np.arange(10000):
        #     smatrix = np.array([[ 0.0,  0,  1,  0,  0,  0,  1,  0,  1,  0,  0,  0],
        #                            [ 0,  1,  0,  1,  0,  1,  0,  0,  0,  0,  1,  1],
        #                            [ 0,  0,  0,  0,  1,  0,  0,  1,  0,  1,  0,  0],
        #                            [ 1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0]])
        #     M = np.array([[ 0.1435,    0.117,     0.0615,    0.0284715, 0.001,     0.0435,    0.001, 0.0085,    0.005,     0.0655,    0.25    ],
        #                        [ 0.248,     0.2425,    0.536,     0.001,     0.0374625, 0.0635,    0.001, 0.021,     0.2,       0.2315,    0.079   ],
        #                        [ 0.348,     0.2335,    0.0745,    0.0034965, 0.935064,  0.035,     0.991513, 0.924,     0.1255,    0.0405,    0.1445  ],
        #                        [ 0.2605,    0.407,     0.328,     0.967032,  0.0264735, 0.858,     0.006487, 0.0465,    0.6695,    0.6625,    0.5265  ]])
        #     #print  smatrix[::,:-1]
        #     [ratio,tarnum,tarsite,tarratio] = core.motif_scan_core(smatrix, M, max_score, score_cutoff, B)
            # print ratio
            # print tarnum
            # print tarsite
            # print tarratio
        for i in np.arange(10000):
            smatrix = np.random.rand(4,1000)
            M = np.random.rand(4,12)
            max_score = 1.2
            score_cutoff = 0.8
            B = [0.28,0.22,0.22,0.28]
            #core.motif_scan_core_2(smatrix, M, max_score, score_cutoff, B)
            core.motif_scan_core(smatrix, M, max_score, score_cutoff, B)

    def test_motif_scan_general(self):
        peak_path = 'Ese14_H3K27ac_LICR_Rep2_peaks_1000.xls'
        peak_format = 'macs'
        #motif_list_path = '/mnt/MAmotif/motifscan_pkg/MotifScan/motif_list/jaspar_2014_vertebrates_motif_list.txt'
        motif_list_path = 'Arnt.txt'
        motif_path = '/mnt/MAmotif/motifscan_pkg/MotifScan/motif/mouse/mouse_jaspar_2014_all_1e-4.txt'
        genome_name = '/home/jiawei/.MotifScan/genome/mm9'
        gene_path = '/home/jiawei/.MotifScan/gene/mm9/refSeq.txt'
        random_times = 5
        peak_length = 1000
        output_dir = 'test_motifscan_general_output'
        region = "genome"
        up = 1000
        down = 1000
        is_enrichment = False
        extract_target_site = True
        motifscan_general.motifscan_general(peak_path,peak_format,motif_path,motif_list_path,genome_name,gene_path,random_times,peak_length,output_dir,region,up,down,is_enrichment,extract_target_site)



    def test_motif_scan(self):
        peak_path = 'Ese14_H3K27ac_LICR_Rep2_peaks.xls'
        p = peak.load_peak(peak_path,'/home/jiawei/.MotifScan/genome/hg19',1000,'macs')
        m = motif.load_motif('/mnt/MAmotif/motifscan_pkg/MotifScan/motif_list/jaspar_2014_vertebrates_motif_list.txt','/mnt/MAmotif/motifscan_pkg/MotifScan/motif/mm9/mm9_jaspar_2014_all_1e-4.txt')
        core.motif_scan(p,m,np.array([0.29,0.21,0.21,0.29]),'./tmp')

    def test_c_motif_scan_core(self):
        #prep
        class MAT(Structure):
            _fields_ = [("n", c_int),
                       ("a_arr",POINTER(c_double)),
                       ("c_arr",POINTER(c_double)),
                       ("g_arr",POINTER(c_double)),
                       ("t_arr",POINTER(c_double))]

        class MOTIF_RES(Structure):
            _fields_ = [("tarnum", c_int),
                       ("ratio",c_double),
                       ("tarsite",POINTER(c_int)),
                       ("tarratio",POINTER(c_double))]

        score_c = CDLL("/usr/local/lib/python2.7/dist-packages/MotifScan/score_c.so")
        score_c.motif_scan_core.restype = POINTER(MOTIF_RES)
        score_c.motif_scan_core.argtypes = [POINTER(MAT),POINTER(MAT),POINTER(c_double*4),c_double,c_double]
        score_c.freeMOTIF_RES.argtypes = [POINTER(MOTIF_RES)]
        #input
        max_score = c_double(8.431914462)
        score_cutoff = c_double(0.797036073733)
        B = (c_double*4)(0.29,0.21,0.21,0.29)

        for i in np.arange(1):
            #np_smatrix = np.random.rand(4,1000)

            np_smatrix = np.array([[ 0,  0,  0,  1,  0,  1],
                      [ 0,  0,  1,  0,  0,  0],
                      [ 1,  0,  0,  0,  1,  0],
                      [ 0,  1,  0,  0,  0,  0]],dtype=np.double)
            a_s = np.ctypeslib.as_ctypes(np_smatrix[0])
            c_s = np.ctypeslib.as_ctypes(np_smatrix[1]) 
            g_s = np.ctypeslib.as_ctypes(np_smatrix[2])  
            t_s = np.ctypeslib.as_ctypes(np_smatrix[3])    

            # np_mat = np.random.rand(4,10)
            np_mat = np.array([[0.1996 , 0.9481 , 0.001 ,  0.001 ,  0.001,   0.001],
                               [0.7984 , 0.001  , 0.997 ,  0.001 ,  0.001,   0.001],
                               [0.001  , 0.0499 , 0.001 ,  0.997 ,  0.001,   0.997],
                               [0.001  , 0.001  , 0.001 ,  0.001 ,  0.997,   0.001]],dtype=np.double)
            a = np.ctypeslib.as_ctypes(np_mat[0])
            c = np.ctypeslib.as_ctypes(np_mat[1])
            g = np.ctypeslib.as_ctypes(np_mat[2])
            t = np.ctypeslib.as_ctypes(np_mat[3])
        
            mat = MAT(6,a,c,g,t)
            smatrix = MAT(6,a_s,c_s,g_s,t_s)
            mat_ptr = score_c.motif_scan_core(byref(smatrix),byref(mat),byref(B),max_score,score_cutoff)

            # ratio = mat_ptr.contents.ratio
            # tarnum = mat_ptr.contents.tarnum
            # b = np.squeeze(np.ctypeslib.as_array(mat_ptr.contents.tarratio,
            #                                 shape=(1,mat_ptr.contents.tarnum)))
            print mat_ptr.contents.ratio
            # print mat_ptr.contents.tarnum
            #tarnum = mat_ptr.contents.tarnum
            #print np.squeeze(np.ctypeslib.as_array(mat_ptr.contents.tarsite,shape=(1,mat_ptr.contents.tarnum)))
            # a = (c_int * 991).from_address(addressof(mat_ptr.contents.tarsite.contents))
            # for i in np.arange(tarnum):
            #     print a[i]
            #
            # print np.squeeze(np.ctypeslib.as_array(mat_ptr.contents.tarratio,shape=(1,mat_ptr.contents.tarnum)))
            # print 'tarsiste:'
            # # for i in np.arange(mat_ptr.contents.tarnum):
            # #     print mat_ptr.contents.tarsite[i]
            # for i in np.arange(10):
            #     print mat_ptr.contents.tarratio[i]

        # print np.ctypeslib.as_array(mat_ptr.contents.tarsite,shape=(1,mat_ptr.contents.tarnum))
        # print np.ctypeslib.as_array(mat_ptr.contents.tarratio,shape=(1,1000-3+1))
    
    def test_c_motif_scan_core_simulation(self):
        #prep
        class MAT(Structure):
            _fields_ = [("n", c_int),
                       ("a_arr",POINTER(c_double)),
                       ("c_arr",POINTER(c_double)),
                       ("g_arr",POINTER(c_double)),
                       ("t_arr",POINTER(c_double))]

        score_c = CDLL("/usr/local/lib/python2.7/dist-packages/score_c.so")
        score_c.motif_scan_core_simulation.restype = c_double
        score_c.motif_scan_core_simulation.argtypes = [POINTER(MAT),POINTER(MAT),POINTER(c_double*4),c_double]
       
        #input
        max_score = c_double(0.8)
        B = (c_double*4)(0.28,0.22,0.22,0.28)

        for i in np.arange(1):
            np_smatrix = np.array([[ 0,  0,  1,  0,  0,  0,  1,  0,  1,  0,  0],
                            [ 0,  1,  0,  1,  0,  1,  0,  0,  0,  0,  1],
                            [ 0,  0,  0,  0,  1,  0,  0,  1,  0,  1,  0],
                            [ 1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0]],dtype=np.double)#np.random.rand(4,10)

            a_s = np.ctypeslib.as_ctypes(np_smatrix[0])
            c_s = np.ctypeslib.as_ctypes(np_smatrix[1]) 
            g_s = np.ctypeslib.as_ctypes(np_smatrix[2])  
            t_s = np.ctypeslib.as_ctypes(np_smatrix[3])    

            np_mat =  np.array([[ 0.1435,    0.117,     0.0615,    0.0284715, 0.001,     0.0435,    0.001, 0.0085,    0.005,     0.0655,    0.25    ],
                      [ 0.248,     0.2425,    0.536,     0.001,     0.0374625, 0.0635,    0.001, 0.021,     0.2,       0.2315,    0.079   ],
                      [ 0.348,     0.2335,    0.0745,    0.0034965, 0.935064,  0.035,     0.991513, 0.924,     0.1255,    0.0405,    0.1445  ],
                      [ 0.2605,    0.407,     0.328,     0.967032,  0.0264735, 0.858,     0.006487, 0.0465,    0.6695,    0.6625,    0.5265  ]]) 
            #np.random.rand(4,10)
            a = np.ctypeslib.as_ctypes(np_mat[0])
            c = np.ctypeslib.as_ctypes(np_mat[1])
            g = np.ctypeslib.as_ctypes(np_mat[2])
            t = np.ctypeslib.as_ctypes(np_mat[3])
        
            mat = MAT(11,a,c,g,t)
            smatrix = MAT(11,a_s,c_s,g_s,t_s)
            print score_c.motif_scan_core_simulation(byref(smatrix),byref(mat),byref(B),max_score)

    def test_simulation_core(self):
        B = np.array([0.28,0.22,0.22,0.28])
        max_score = 0.8
        smatrix = np.array([[ 0,  0,  1,  0,  0,  0,  1,  0,  1,  0,  0],
                            [ 0,  1,  0,  1,  0,  1,  0,  0,  0,  0,  1],
                            [ 0,  0,  0,  0,  1,  0,  0,  1,  0,  1,  0],
                            [ 1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0]],dtype=np.double)
        M = np.array([[ 0.1435,    0.117,     0.0615,    0.0284715, 0.001,     0.0435,    0.001, 0.0085,    0.005,     0.0655,    0.25    ],
                      [ 0.248,     0.2425,    0.536,     0.001,     0.0374625, 0.0635,    0.001, 0.021,     0.2,       0.2315,    0.079   ],
                      [ 0.348,     0.2335,    0.0745,    0.0034965, 0.935064,  0.035,     0.991513, 0.924,     0.1255,    0.0405,    0.1445  ],
                      [ 0.2605,    0.407,     0.328,     0.967032,  0.0264735, 0.858,     0.006487, 0.0465,    0.6695,    0.6625,    0.5265  ]])
        print motif.simulation_core(smatrix,M,max_score,B)

    def test_c_copyMAT(self):
        class MAT(Structure):
            _fields_ = [("n", c_int),
                       ("a_arr",POINTER(c_double)),
                       ("c_arr",POINTER(c_double)),
                       ("g_arr",POINTER(c_double)),
                       ("t_arr",POINTER(c_double))]
        score_c = CDLL("../lib/score.so")
        score_c.logMAT.restype = POINTER(MAT)
        score_c.logMAT.argtypes = [POINTER(MAT)]
        a = (c_double*3)(0.2,0.2,0.2)
        c = (c_double*3)(0.1,0.3,0.2)
        g = (c_double*3)(0.4,0.2,0.5)
        t = (c_double*3)(0.3,0.3,0.1)
        a_ptr = cast(pointer(a),POINTER(c_double))
        c_ptr = cast(pointer(c),POINTER(c_double))
        g_ptr = cast(pointer(g),POINTER(c_double))
        t_ptr = cast(pointer(t),POINTER(c_double))
        mat = MAT(3,a_ptr,c_ptr,g_ptr,t_ptr)
        m_ptr = score_c.logMAT(pointer(mat)) 

    def test_peak2peak_enrichment(self):
        res1 = pd.read_csv('motifscan_output_Ese14_H3K27ac_LICR_Rep2_peaks/motif_enrichment.csv',index_col=0)
        res2 = pd.read_csv('motifscan_output_peakSeq.optimal.wgEncodeSydhTfbsK562bZnf263UcdAlnRep0_vs_wgEncodeSydhTfbsK562bInputUcdAlnRep1_new_summit/motif_enrichment.csv',index_col=0)

        motif_table = core.target_enrichment_peak2peak(res1, res2)
        motif_table.to_csv("peak2peak_enrichment.csv")

    def test_deduplicate(self):
        ts = [1, 5, 6, 7, 18,30]
        tr = [0.5, 0.3, 0.3, 0.2, 0.7, 0.3]
        motif_len = 5
        ts, tr = core.deduplicate_target_site(ts, tr, motif_len)
        print ts, tr

    def test_calculate_ES(self):
        peak_table = pd.read_pickle('test_peak_table.pkl')
        target_yes, target_no = core.split_manorm_peaks_motif(peak_table, 'TAL1::GATA1')
        target_yes_idx =  target_yes.index
        L = peak_table['value'].copy().squeeze()
        L.sort()
     
        ES, hit_cum = core.calculate_ES(L,target_yes_idx)
        core.ES_plot(hit_cum,L,target_yes_idx,'psea_test.png')

    def test_peak_set_enrichemnt_analysis(self):
        # peak_table = pd.read_pickle('test_peak_table.pkl')
        p1 = pd.read_pickle('MAmotif_test_input/h1hesc_peak_result.pkl')
        r1 = pd.read_pickle('MAmotif_test_input/h1hesc_rnd_result.pkl')
        e1 = pd.read_pickle('MAmotif_test_input/h1hesc_enrich_result.pkl')
        p2 = pd.read_pickle('MAmotif_test_input/k562_peak_result.pkl')
        r2 = pd.read_pickle('MAmotif_test_input/k562_rnd_result.pkl')
        e2 = pd.read_pickle('MAmotif_test_input/k562_enrich_result.pkl')
        p3 = pd.read_pickle('MAmotif_test_input/common_peak_result.pkl')
        r3 = pd.read_pickle('MAmotif_test_input/common_rnd_result.pkl')
        e3 = pd.read_pickle('MAmotif_test_input/common_enrich_result.pkl')

        [p,r,e] = motifscan_general.merge_two_results(p1,p2,r1,r2,e1)
        [p,r,e] = motifscan_general.merge_two_results(p,p3,r,r3,e3)
        core.peak_set_enrichment_analysis(p,'test_fc_plot/')

    def test_peak_result_test(self):
        p1 = pd.read_pickle('MAmotif_test_input/h1hesc_peak_result.pkl')
        r1 = pd.read_pickle('MAmotif_test_input/h1hesc_rnd_result.pkl')
        e1 = pd.read_pickle('MAmotif_test_input/h1hesc_enrich_result.pkl')
        p2 = pd.read_pickle('MAmotif_test_input/k562_peak_result.pkl')
        r2 = pd.read_pickle('MAmotif_test_input/k562_rnd_result.pkl')
        e2 = pd.read_pickle('MAmotif_test_input/k562_enrich_result.pkl')
        p3 = pd.read_pickle('MAmotif_test_input/common_peak_result.pkl')
        r3 = pd.read_pickle('MAmotif_test_input/common_rnd_result.pkl')
        e3 = pd.read_pickle('MAmotif_test_input/common_enrich_result.pkl')
        [p,r,e] = motifscan_general.merge_two_results(p1,p3,r1,r3,e1)
        core.peak_result_test(p2)

    def test_split_MAmotif(self):
        self.peak_table = pd.read_pickle('test_peak_table.pkl')
        self.motif_table = pd.read_pickle('test_motif_table.pkl')
        core.split_manorm_peaks_motif(self.peak_table, self.motif_table)

    def test_h_cluster(self):
        p1 = pd.read_pickle('MAmotif_test_input/h1hesc_peak_result.pkl')
        r1 = pd.read_pickle('MAmotif_test_input/h1hesc_rnd_result.pkl')
        e1 = pd.read_pickle('MAmotif_test_input/h1hesc_enrich_result.pkl')
        p2 = pd.read_pickle('MAmotif_test_input/k562_peak_result.pkl')
        r2 = pd.read_pickle('MAmotif_test_input/k562_rnd_result.pkl')
        e2 = pd.read_pickle('MAmotif_test_input/k562_enrich_result.pkl')
        [p,r,e] = motifscan_general.merge_two_results(p1,p2,r1,r2,e1)
        core.h_cluster(p,e,'./test_fc_plot',1.3,0.001)

    def test_fc_plot(self):
        p1 = pd.read_pickle('MAmotif_test_input/h1hesc_peak_result.pkl')
        r1 = pd.read_pickle('MAmotif_test_input/h1hesc_rnd_result.pkl')
        e1 = pd.read_pickle('MAmotif_test_input/h1hesc_enrich_result.pkl')
        p2 = pd.read_pickle('MAmotif_test_input/k562_peak_result.pkl')
        r2 = pd.read_pickle('MAmotif_test_input/k562_rnd_result.pkl')
        e2 = pd.read_pickle('MAmotif_test_input/k562_enrich_result.pkl')
        p3 = pd.read_pickle('MAmotif_test_input/common_peak_result.pkl')
        r3 = pd.read_pickle('MAmotif_test_input/common_rnd_result.pkl')
        e3 = pd.read_pickle('MAmotif_test_input/common_enrich_result.pkl')

        [p,r,e] = motifscan_general.merge_two_results(p1,p2,r1,r2,e1)
        [p,r,e] = motifscan_general.merge_two_results(p,p3,r,r3,e3)

        core.fc_plot(p,r,'./test_fc_plot')

    def test_target_site_distribution(self):
        peak_table = pd.read_pickle("test_peak_table.pkl")
        motif_table = pd.read_pickle("test_motif_table.pkl")
        core.target_site_distribution(peak_table, motif_table, "test_motif_target_plot", 0)

    def test_c_sliding_score(self):
        # s = (c_double*5)(1,2,3,4,5)
        # w = (c_double*3)(1,2,3)
        s = np.array([1,2,3,4,5]).ctypes.data_as(POINTER(c_double))
        w = np.array([1,2,3]).ctypes.data_as(POINTER(c_double))

        score_c = CDLL("../lib/score.so")
        score_c.sliding_score.restype = c_double
        print score_c.sliding_score(w,3,s,5)

    def test_ctype_structure(self):
        class MAT(Structure):
            _fields_ = [("n", c_int),
                       ("a",POINTER(c_double)),
                       ("c",POINTER(c_double)),
                       ("g",POINTER(c_double)),
                       ("t",POINTER(c_double))]
        a = (c_double*4)(1,2,3,4)
        c = (c_double*4)(1,2,3,4)
        g = (c_double*4)(1,2,3,4)
        t = (c_double*4)(1,2,3,4)


        mat = MAT(4,a,c,g,t)
        win = MAT(4,a,c,g,t)
        score_c = CDLL("../lib/score.so")
        score_c.sliding_score_mat.restype = c_double
        print score_c.sliding_score_mat(mat,win)

    def test_peak_test(self):
        common_motif_dir = '../MAmotif_H3K4me3_H1hesc_VS_K562/motifscan_on_common_peaks'
        a_motif_dir = '../MAmotif_H3K4me3_H1hesc_VS_K562/motifscan_on_H1hesc_H3k4me3_Broad_Rep1_P100_unique'

        print 'Loading motifscan results...'
        common_peak_result = pd.read_pickle("%s/peak_result.pkl"%common_motif_dir)
        common_rnd_result = pd.read_pickle("%s/rnd_result.pkl"%common_motif_dir)
        common_enrich_result = pd.read_pickle("%s/enrich_result.pkl"%common_motif_dir)
        a_peak_result = pd.read_pickle("%s/peak_result.pkl"%a_motif_dir)
        a_rnd_result = pd.read_pickle("%s/rnd_result.pkl"%a_motif_dir)
        
        print 'Merging peaks...'
        [positive_peak_result,positive_rnd_result,positive_enrich_result] = \
    motifscan_general.merge_two_results(common_peak_result,a_peak_result,common_rnd_result,a_rnd_result,common_enrich_result)

        core.peak_result_test(positive_peak_result, positive_rnd_result,plot_output_dir='enrichment_plot')

    def test_simulation(self):
        raw_motif_path = '../example/jaspar_vertebrata.txt'
        #raw_motif_path = '/home/jiawei/yijing_motif/motif/CArG_motif_matrix_t.txt'
        genome_db_path = '/home/jiawei/.MotifScan/genome/tair10'
        sample_number = 100000
        background = pd.read_pickle('%s/background'%genome_db_path)['background']
        chromosome_size = pd.read_pickle('%s/chromosome_size'%genome_db_path)
        motif_table = motif.load_motif_matrix(raw_motif_path)
        motif_table = motif.compute_max_score(motif_table,background)
        motif_table = motif.simulation(motif_table,sample_number,genome_db_path,chromosome_size,background)
        print motif_table["score_cutoff_3"]
    
    def test_compute_genome_background(self):
        #genome_path = "/opt/genome_collection/ce10/ce10.fa"
        genome_path = "/opt/genome_collection/sacCer3/sacCer3.fa"
        genome.compute_genome_background(genome_path)

    def test_write_motif_table(self):
        motif_table = pd.read_pickle("/home/jiawei/.MotifScan/motif/hg19/jaspar_all_motif_hg19_20150209134637.pkl")
        fo = open('motif_db.txt','w')
        motif.write_motif_table(motif_table,fo)

    def test_read_motif_table(self):
        fi = open('../motif/hg19/hg19_jaspar_2014_all.txt','r')
        motif_table = motif.read_motif_table(fi)
        print motif_table.iloc[0]['matrix']
    
    def test_load_motif(self):
        #motif_db_path = '/home/jiawei/MAmotif/motifscan_pkg/MotifScan/tests/test_motif_db'
        motif_path = '/mnt/MAmotif/motifscan_pkg/MotifScan/motif/hg19/hg19_jaspar_2014_all_1e-4.txt'
        motif_list_path = ''
        motif_table = motif.load_motif(motif_list_path, motif_path)
    
    def test_load_motif_matrix(self):
        matrix_file = '/home/jiawei/MAmotif/motifscan_pkg/MotifScan/example/jaspar_all_motif.txt'
        motif.load_motif_matrix(matrix_file)

def suite():
    suite = unittest.TestSuite()
    #suite.addTest(testMotifScanCore("test_sliding_score"))
    #suite.addTest(testMotifScanCore("test_sliding_score_mat"))
    #suite.addTest(testMotifScanCore("test_motif_scan_core"))

    #suite.addTest(testMotifScanCore("test_c_motif_scan_core")) 

    #suite.addTest(testMotifScanCore("test_c_motif_scan_core_simulation"))
    #suite.addTest(testMotifScanCore("test_simulation_core"))

    #suite.addTest(testMotifScanCore("test_motif_scan"))
    #suite.addTest(testMotifScanCore("test_peak2peak_enrichment"))
    #suite.addTest(testMotifScanCore("test_motif_scan_general"))
    #suite.addTest(testMotifScanCore("test_c_copyMAT"))
    #suite.addTest(testMotifScanCore("test_c_sliding_score_mat"))
    #suite.addTest(testMotifScanCore("test_c_testMAT"))
    #suite.addTest(testMotifScanCore("test_c_testMOTIF_RES"))
    #suite.addTest(testMotifScanCore("test_peak_set_enrichemnt_analysis"))
    #suite.addTest(testMotifScanCore("test_calculate_ES"))
    #suite.addTest(testMotifScanCore("test_peak_result_test"))
    #suite.addTest(testMotifScanCore("test_split_MAmotif"))
    #suite.addTest(testMotifScanCore("test_tarsite_fc_plot"))
    #suite.addTest(testMotifScanCore("test_fc_plot"))
    suite.addTest(testMotifScanCore("test_target_site_distribution"))
    #suite.addTest(testMotifScanCore("test_c_sliding_score"))
    #suite.addTest(testMotifScanCore("test_ctype_structure"))
    #suite.addTest(testMotifScanCore("test_peak_test"))
    #suite.addTest(testMotifScanCore("test_h_cluster"))
    #suite.addTest(testMotifScanCore("test_simulation"))
    #suite.addTest(testMotifScanCore("test_load_motif"))
    #suite.addTest(testMotifScanCore("test_compute_genome_background"))
    #suite.addTest(testMotifScanCore("test_write_motif_table"))
    #suite.addTest(testMotifScanCore("test_read_motif_table"))
    #suite.addTest(testMotifScanCore("test_load_motif"))
    #suite.addTest(testMotifScanCore("test_deduplicate"))
    #suite.addTest(testMotifScanCore("test_load_motif_matrix"))
    return suite

if __name__ == '__main__':
    unittest.main(defaultTest='suite')