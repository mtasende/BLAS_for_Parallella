% 
% --- BLIS library info -------------------------------------
% 
% version string               0.1.6
% 
% --- BLIS config header ---
% 
% integer type size (bits)     32
% # of floating-point types    4
% maximum type size            16
% 
% maximum number of threads    1
% 
% SIMD alignment (bytes)       16
% 
% stack memory allocation        
%   address alignment (bytes)  16
% 
% dynamic memory allocation      
%   address alignment          16
%   stride alignment           64
% 
% contiguous memory allocation   
%   # of mc x kc blocks        1
%   # of kc x nc blocks        1
%   # of mc x nc blocks        0
%   block address alignment    4096
%   max preload byte offset    128
%   actual pool sizes (bytes)    
%     for mc x kc blocks of A  21256320
%     for kc x nc panels of B  21256320
%     for mc x nc panels of C  128
% 
% BLAS compatibility layer       
%   enabled?                   1
%   integer type size (bits)   32
% 
% floating-point types           s       d       c       z 
%   sizes (bytes)                4       8       8      16
% 
% 
% --- BLIS default implementations ---
% 
% level-3 implementations        s       d       c       z
%   gemm                    native  native  native  native
%   hemm                    native  native  native  native
%   herk                    native  native  native  native
%   her2k                   native  native  native  native
%   symm                    native  native  native  native
%   syrk                    native  native  native  native
%   syr2k                   native  native  native  native
%   trmm                    native  native  native  native
%   trmm3                   native  native  native  native
%   trsm                    native  native  native  native
% 
% --- BLIS native implementation info ---
% 
%                                                c       z 
% complex implementation                    native  native
% 
% level-3 blocksizes             s       d       c       z 
%   mc                         768     768     128     128
%   kc                        4096    2048     256     256
%   nc                         768     768    2048    2048
% 
%   mc maximum                 768     768     128     128
%   kc maximum                4096    2048     256     256
%   nc maximum                 768     768    2048    2048
% 
%   mr                         192     192      16      16
%   nr                         256     256      16      16
% 
%   mr packdim                 192     192      16      16
%   nr packdim                 256     256      16      16
% 
% micro-kernel types             s       d       c       z
%   gemm                   optimzd optimzd optimzd optimzd
%   gemmtrsm_l             optimzd optimzd optimzd optimzd
%   gemmtrsm_u             optimzd optimzd optimzd optimzd
%   trsm_l                 optimzd optimzd optimzd optimzd
%   trsm_u                 optimzd optimzd optimzd optimzd
% 
% 
% --- BLIS induced implementation info ---
% 
%                                c       z 
% complex implementation       3mh     3mh
% 
% level-3 blocksizes             c       z 
%   mc                         768     768
%   kc                        4096    2048
%   nc                         768     768
% 
%   mc maximum                 768     768
%   kc maximum                4096    2048
%   nc maximum                 768     768
% 
%   mr                         192     192
%   nr                         256     256
% 
%   mr packdim                 192     192
%   nr packdim                 256     256
% 
% micro-kernel types             c       z
%   gemm                   virtual virtual
%   gemmtrsm_l             notappl notappl
%   gemmtrsm_u             notappl notappl
%   trsm_l                 notappl notappl
%   trsm_u                 notappl notappl
% 
%                                c       z 
% complex implementation       3m3     3m3
% 
% level-3 blocksizes             c       z 
%   mc                         768     768
%   kc                        4096    2048
%   nc                         256     256
% 
%   mc maximum                 768     768
%   kc maximum                4096    2048
%   nc maximum                 256     256
% 
%   mr                         192     192
%   nr                         256     256
% 
%   mr packdim                 192     192
%   nr packdim                 256     256
% 
% micro-kernel types             c       z
%   gemm                   virtual virtual
%   gemmtrsm_l             notappl notappl
%   gemmtrsm_u             notappl notappl
%   trsm_l                 notappl notappl
%   trsm_u                 notappl notappl
% 
%                                c       z 
% complex implementation       3m2     3m2
% 
% level-3 blocksizes             c       z 
%   mc                         192     192
%   kc                        4096    2048
%   nc                         256     256
% 
%   mc maximum                 192     192
%   kc maximum                4096    2048
%   nc maximum                 256     256
% 
%   mr                         192     192
%   nr                         256     256
% 
%   mr packdim                 192     192
%   nr packdim                 256     256
% 
% micro-kernel types             c       z
%   gemm                   virtual virtual
%   gemmtrsm_l             notappl notappl
%   gemmtrsm_u             notappl notappl
%   trsm_l                 notappl notappl
%   trsm_u                 notappl notappl
% 
%                                c       z 
% complex implementation       3m1     3m1
% 
% level-3 blocksizes             c       z 
%   mc                         768     768
%   kc                        1365     682
%   nc                         768     768
% 
%   mc maximum                 768     768
%   kc maximum                1365     682
%   nc maximum                 768     768
% 
%   mr                         192     192
%   nr                         256     256
% 
%   mr packdim                 192     192
%   nr packdim                 256     256
% 
% micro-kernel types             c       z
%   gemm                   virtual virtual
%   gemmtrsm_l             virtual virtual
%   gemmtrsm_u             virtual virtual
%   trsm_l                 virtual virtual
%   trsm_u                 virtual virtual
% 
%                                c       z 
% complex implementation       4mh     4mh
% 
% level-3 blocksizes             c       z 
%   mc                         768     768
%   kc                        4096    2048
%   nc                         768     768
% 
%   mc maximum                 768     768
%   kc maximum                4096    2048
%   nc maximum                 768     768
% 
%   mr                         192     192
%   nr                         256     256
% 
%   mr packdim                 192     192
%   nr packdim                 256     256
% 
% micro-kernel types             c       z
%   gemm                   virtual virtual
%   gemmtrsm_l             notappl notappl
%   gemmtrsm_u             notappl notappl
%   trsm_l                 notappl notappl
%   trsm_u                 notappl notappl
% 
%                                c       z 
% complex implementation      4m1b    4m1b
% 
% level-3 blocksizes             c       z 
%   mc                         384     384
%   kc                        4096    2048
%   nc                         256     256
% 
%   mc maximum                 384     384
%   kc maximum                4096    2048
%   nc maximum                 256     256
% 
%   mr                         192     192
%   nr                         256     256
% 
%   mr packdim                 192     192
%   nr packdim                 256     256
% 
% micro-kernel types             c       z
%   gemm                   virtual virtual
%   gemmtrsm_l             notappl notappl
%   gemmtrsm_u             notappl notappl
%   trsm_l                 notappl notappl
%   trsm_u                 notappl notappl
% 
%                                c       z 
% complex implementation      4m1a    4m1a
% 
% level-3 blocksizes             c       z 
%   mc                         768     768
%   kc                        2048    1024
%   nc                         768     768
% 
%   mc maximum                 768     768
%   kc maximum                2048    1024
%   nc maximum                 768     768
% 
%   mr                         192     192
%   nr                         256     256
% 
%   mr packdim                 192     192
%   nr packdim                 256     256
% 
% micro-kernel types             c       z
%   gemm                   virtual virtual
%   gemmtrsm_l             virtual virtual
%   gemmtrsm_u             virtual virtual
%   trsm_l                 virtual virtual
%   trsm_u                 virtual virtual
% 
% 
% --- BLIS misc. other info ---
% 
% micro-panel alignment (bytes)  s       d       c       z 
%   A (left matrix)              4       8       8      16
%   B (right matrix)             4       8       8      16
% 
% level-2 cache blocksizes       s       d       c       z 
%   m dimension               1000    1000    1000    1000
%   n dimension               1000    1000    1000    1000
% 
% level-1f fusing factors        s       d       c       z 
%   default                      8       4       4       2
%   axpyf                        8       4       4       2
%   dotxf                        8       4       4       2
%   dotxaxpyf                    8       4       4       2
% 

% 
% --- BLIS test suite parameters ----------------------------
% 
% num repeats per experiment   1
% num matrix storage schemes   1
% storage[ matrix ]            c
% num vector storage schemes   1
% storage[ vector ]            c
% mix all storage schemes?     1
% general stride spacing       2048
% num datatypes                1
% datatype[0]                  0 (s)
% problem size: first to test  768
% problem size: max to test    768
% problem size increment       1
% test induced complex           
%   3mh?                       1
%   3m3?                       1
%   3m2?                       1
%   3m1?                       1
%   4mh?                       1
%   4m1b (4mb)?                1
%   4m1a (4m1)?                1
% test native complex?         1
% error-checking level         1
% reaction to failure          i
% output in matlab format?     1
% output to stdout AND files?  1
% 

% --- gemm ---
% 
% test gemm seq front-end?    1
% gemm m n k                  4096 4096 4096
% gemm operand params         ??
% 

% blis_<dt><op>_<params>_<stor>      m     n     k   gflops   resid      result
blis_sgemm_nn_ccc         (   1, 1:5 ) = [  4096  4096  4096    2.373  4.52e-07 ]; % PASS
blis_sgemm_nc_ccc         (   1, 1:5 ) = [  4096  4096  4096    2.372  4.79e-07 ]; % PASS
blis_sgemm_nt_ccc         (   1, 1:5 ) = [  4096  4096  4096    2.447  4.77e-07 ]; % PASS
blis_sgemm_nh_ccc         (   1, 1:5 ) = [  4096  4096  4096    2.447  4.65e-07 ]; % PASS
blis_sgemm_cn_ccc         (   1, 1:5 ) = [  4096  4096  4096    2.368  4.69e-07 ]; % PASS
blis_sgemm_cc_ccc         (   1, 1:5 ) = [  4096  4096  4096    2.368  4.75e-07 ]; % PASS
blis_sgemm_ct_ccc         (   1, 1:5 ) = [  4096  4096  4096    2.446  4.67e-07 ]; % PASS
blis_sgemm_ch_ccc         (   1, 1:5 ) = [  4096  4096  4096    2.445  4.59e-07 ]; % PASS
blis_sgemm_tn_ccc         (   1, 1:5 ) = [  4096  4096  4096    2.019  4.50e-07 ]; % PASS
blis_sgemm_tc_ccc         (   1, 1:5 ) = [  4096  4096  4096    2.022  4.64e-07 ]; % PASS
blis_sgemm_tt_ccc         (   1, 1:5 ) = [  4096  4096  4096    2.074  4.55e-07 ]; % PASS
blis_sgemm_th_ccc         (   1, 1:5 ) = [  4096  4096  4096    2.075  4.89e-07 ]; % PASS
blis_sgemm_hn_ccc         (   1, 1:5 ) = [  4096  4096  4096    2.023  4.67e-07 ]; % PASS
blis_sgemm_hc_ccc         (   1, 1:5 ) = [  4096  4096  4096    2.015  4.69e-07 ]; % PASS
blis_sgemm_ht_ccc         (   1, 1:5 ) = [  4096  4096  4096    2.077  4.69e-07 ]; % PASS
blis_sgemm_hh_ccc         (   1, 1:5 ) = [  4096  4096  4096    2.070  4.63e-07 ]; % PASS

