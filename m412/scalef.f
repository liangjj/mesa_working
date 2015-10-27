*deck @(#)scalef.f	5.1  11/6/94
      subroutine scalef(ex,lensto,atomno,shlpnt)
c***begin prologue     scalef
c***date written       850601  yymmdd
c***revision date      yymmdd  yymmdd
c***keywords           minimum basis, scale factors, basis
c***author             martin, richard (lanl)
c***source             @(#)scalef.f	5.1   11/6/94
c***purpose            returns scale factor for a specific atomic orbital/atom.
c***description
c     call scalef(ex,lensto,atomno,shlpnt)
c       ex      the exponent vector(lensto).  the input exponents
c               are scaled upon return by the scalefactor**2.
c       lensto  the expansion length.
c       atomno  the atomic number.
c       shlpnt  the orbital type.
c***references         (none)
c***routines called    defsc(m412)
c***end prologue       scalef
      implicit integer(a-z)
      real*8 ex(lensto),defsc,zero,scfac
      real*8 scale(110,19),scal1s(110),scal2s(110),scal2p(110),
     $          scal3s(110), scal3p(110), scal4s(110), scal4p(110),
     $          scal3d(110), scal5s(110), scal5p(110), scal4d(110),
     $          scal6s(110), scal6p(110), scal5d(110), scal4f(110),
     $          scal7s(110), scal7p(110), scal6d(110), scal5f(110)
      equivalence (scal1s(1),scale(1,1)), (scal2s(1),scale(1,2)),
     $            (scal2p(1),scale(1,3)), (scal3s(1),scale(1,4)),
     $            (scal3p(1),scale(1,5)), (scal4s(1),scale(1,6)),
     $            (scal4p(1),scale(1,7)), (scal3d(1),scale(1,8)),
     $            (scal5s(1),scale(1,9)), (scal5p(1),scale(1,10)),
     $            (scal4d(1),scale(1,11)), (scal6s(1),scale(1,12)),
     $            (scal6p(1),scale(1,13)), (scal5d(1),scale(1,14)),
     $            (scal4f(1),scale(1,15)), (scal7s(1),scale(1,16)),
     $            (scal7p(1),scale(1,17)), (scal6d(1),scale(1,18)),
     $            (scal5f(1),scale(1,19))
      parameter (zero=0.0d00)
c
      data scal1s/
     $          1.24d0,  1.69d0,  2.69d0,  3.68d0,  4.68d0,  5.67d0,
     $          6.67d0,  7.66d0,  8.65d0,  9.64d0, 10.61d0, 11.59d0,
     $         12.56d0, 13.53d0, 14.50d0, 15.47d0, 16.43d0, 17.40d0,
     $    18.4895,19.4730,20.4566,21.4409,22.4256,
     $    23.4138,24.3957,25.3810,26.3668,29.3245,28.3386,29.3245,
     $    30.3094,31.2937,32.2783,33.2622,34.2471,
     $    35.2316,36.2078,37.1911,38.1756,39.1590,40.1423,41.1256,
     $    42.1090,43.0923,44.0756,45.0589,46.0423,47.0256,48.0097,
     $    48.9920,49.9744,50.9568,51.9391,52.9215,56*0./
      data scal2s/0.0d0,0.0d0,
     $    0.65d0,  0.98d0,  1.30d0,  1.63d0,
     $    1.95d0,  2.28d0,  2.43d0,  2.70d0,
     $    3.48d0,  3.90d0,  4.36d0,  4.83d0,
     $    5.31d0,  5.79d0,  6.26d0,  6.74d0,
     $    6.5031,6.8882,7.2868,7.6883,8.0907,8.4919,8.8969,
     $    9.2995,9.7025,10.1063,10.5099,10.9140,11.2995,11.6824,12.0635,
     $    12.4442,12.8217,13.1990,13.5784,13.9509,14.3111,14.6869,
     $    15.0626,15.4384,15.8141,16.1899,16.5656,16.9414,17.3171,
     $    17.6929,18.0618,18.4297,18.7977,19.1656,19.5335,19.9015,56*0./
      data scal2p/2*0.0d0,
     $    0.65d0,  0.98d0,  1.30d0,  1.63d0,
     $    1.95d0,  2.28d0,  2.43d0,  2.70d0,
     $    3.48d0,  3.90d0,  4.36d0,  4.83d0,
     $    5.31d0,  5.79d0,  6.26d0,  6.74d0,
     $    7.5136,8.0207,8.5273,9.0324,9.5364,10.0376,10.5420,11.0444,
     $    11.5462,12.0476,12.5485,13.0490,13.5454,14.0411,14.5368,
     $    15.0326,15.5282,16.0235,16.5194,17.0152,17.5016,17.9964,
     $    18.4911,18.9859,19.4704,19.9754,20.4702,20.9650,21.4597,
     $    21.9545,22.4490,22.9427,23.4363,23.9300,24.4237,24.9173,56*0./
      data scal3s/10*0.0d0,
     $          0.73d0,  0.95d0,  1.17d0,  1.38d0,
     $          1.60d0,  1.82d0,  2.03d0,  2.20d0,
     $    2.8933,3.2005,3.4466,3.6777,3.9031,4.1226,4.3393,4.5587,
     $    4.7741,4.9870,5.1981,5.4064,5.6654,5.9299,6.1985,6.4678,
     $    6.7395,7.0109,7.2809,7.5546,7.8505,8.1205,8.3905,8.6605,
     $    8.9304,9.2004,9.4795,9.7404,10.0104,10.2804,10.5436,
     $    10.8066,11.0697,11.3327,11.5958,11.8588,56*0./
      data scal3p/10*0.0d0,
     $          0.73d0,  0.95d0,  1.17d0,  1.38d0,
     $          1.60d0,  1.82d0,  2.03d0,  2.20d0,
     $    2.8933,3.2005,
     $    3.1354,3.3679,3.5950,3.8220,4.0364,4.2593,4.4782,4.6950,
     $    4.9102,5.1231,5.4012,5.6712,5.9499,6.2350,6.5236,6.8114,
     $    7.1011,7.3892,7.6975,7.9485,8.2052,8.4912,8.7947,9.0737,
     $    9.3848,9.6732,9.9362,10.2305,10.5069,10.7844,11.0613,
     $    11.3363,11.6138,11.8892,56*0./
      data scal4s/18*0.0d0,
     $    0.8738,1.0995,1.1581,1.2042,1.2453,1.2833,1.3208,1.3585,
     $    1.3941,1.4277,1.4606,1.4913,1.7667,2.0109,2.2360,2.4394,
     $    2.6382,2.8289,3.0970,3.3611,3.5659,3.7254,3.8207,4.0241,
     $    4.2996,4.4140,4.6454,4.7465,4.9662,5.2173,5.4403,5.6645,
     $    5.8859,6.1021,6.3243,6.5432,56*0./
      data scal3d/20*0.0d0,
     $    2.3733,2.7138,2.9943,3.2522,3.5094,3.7266,3.9518,4.1765,
     $    4.4002,4.6261,5.0311,5.4171,5.7928,6.1590,6.5197,6.8753,
     $    7.2264,7.5754,8.4657,8.5223,8.7490,9.0761,9.4510,
     $    9.7863,10.1350,10.4837,10.8466,11.2023,11.5594,
     $    11.9139,12.2666,12.6131,12.9669,13.3156,56*0./
      data scal4p/18*0.0d0,
     $    0.8738,1.0995,1.1581,1.2042,1.2453,1.2833,1.3208,1.3585,
     $    1.3941,1.4277,1.4606,1.4913,
     $    1.5554,1.6951,1.8623,2.0718,2.2570,2.4423,2.7202,2.9830,
     $    3.1864,3.3650,3.5211,3.7442,3.9528,4.1087,4.2849,4.4308,
     $    4.6406,4.8528,5.0922,5.3163,5.5453,5.7805,
     $    6.0074,6.2393,56*0./
      data scal5s/36*0.0d0,
     $    0.9969,1.2141,1.2512,1.2891,1.1842,1.2212,
     $    1.4453,1.2969,1.3279,1.5675,1.3511,1.6384,1.9023,
     $    2.1257,2.3222,2.5076,2.6807,2.8436,56*0./
      data scal4d/38*0.0d0,
     $    3.9896,3.2679,2.8094,2.8481,3.2205,3.2032,3.3606,3.4044,
     $    3.6907,3.9692,4.2354,4.4925,4.7436,4.9900,5.2335,5.4733,56*0./
      data scal5p/36*0.0d0,
     $    0.9969,1.2141,1.2512,1.2891,1.1842,1.2212,
     $    1.4453,1.2969,1.3279,1.5675,1.3511,1.6384,
     $    1.6940,1.8204,1.9989,2.1617,2.3223,2.4849,56*0./
      data scal6s/54*0.0d0,
     $    1.0605,1.2625,1.5520,1.7994,1.2911,1.5511,1.5659,1.3353,
     $    1.3536,1.3691,1.3834,1.3906,1.4065,1.4127,1.4307,1.4322,
     $    1.4674,1.5274,1.5875,1.6424,1.6860,1.7205,1.7611,1.7919,
     $    1.8230,1.8589,2.1366,2.3500,2.5400,2.7218,2.8833,3.0540,
     $    24*0.0d0/
      data scal6p/54*0.0d0,
     $    1.0605,1.2625,1.5520,1.7994,1.2911,1.5511,1.5659,1.3353,
     $    1.3536,1.3691,1.3834,1.3906,1.4065,1.4127,1.4307,1.4322,
     $    1.4674,1.5274,1.5875,1.6424,1.6860,1.7205,1.7611,1.7919,
     $    1.8230,1.8589,2.0423,2.0655,2.2233,2.3701,2.5272,2.6793,
     $    24*0.0d0/
      data scal5d/56*0.0d0,
     $    14*4.0d0,4.0226,3.3239,3.2736,3.3484,3.4766,3.5994,3.7392,
     $    3.8815,4.0253,4.1712,4.4050,4.6304,4.8488,5.0608,5.2678,
     $    5.4706,24*0.0d0/
      data scal4f/56*0.d0,
     $    3.5230,4.7000,5.2752,5.5665,5.7835,5.8829,6.0800,6.2534,
     $    6.4662,6.6340,6.8674,6.9946,7.1585,7.3580,7.7328,8.0524,
     $    8.3676,8.6777,8.9812,9.2882,9.5862,9.8765,10.1624,10.4402,
     $    10.7169,10.9922,11.2673,11.5396,11.8101,12.0828,24*0.0d0/
      data scal5f/89*0.0d0,21*1.0d+00/
      data scal6d/89*0.0d0,21*1.0d+00/
      data scal7s/89*0.0d0,21*1.0d+00/
      data scal7p/89*0.0d0,21*1.0d+00/
      save scal1s,scal2s,scal2p,scal3s,scal3p,scal4s,scal3d,scal4p
      save scal5s,scal4d,scal5p,scal6s,scal6p,scal5d,scal4f,scal5f
      save scal6d,scal7s,scal7p
c
c     this routine returns a scale factor for the atomic number atomno,
c     and orbital type shlpnt.
c     shlpnt maps a principal quantum number and angular momentum
c     into an integer.  the values are:
c      1 1s
c      2 2s
c      3 2p
c      4 3s
c      5 3p
c      6 4s
c      7 4p
c      8 3d
c      9 5s
c     10 5p
c     11 4d
c     12 6s
c     13 6p
c     14 5d
c     15 4f
c     16 7s
c     17 7p
c     18 6d
c     19 5f
c
c
c
      scfac=scale(atomno,shlpnt)
      if(scfac.eq.zero) scfac=defsc(atomno,shlpnt)
      do 10 j=1,lensto
         ex(j)=ex(j)*scfac*scfac
   10 continue
c
c
      return
      end