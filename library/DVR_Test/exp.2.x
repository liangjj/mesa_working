 mesa(1.2);5/30/91;unix.                                                        
     (c) 1990, the university of california.                                    
     P.W. Saxe, B.H. Lengsfield iii, R.L. Martin, M. Page and B. Schneider.     

     23-jun-1900 16:23:27    
     National Science Foundation;Dec Alpha 600MHz(bohr);B. Schneider                 
 main files/memory:
     inp:   inp.exp                       out:    mesa.out                 
     chk:   mesa.chk                      dat:    ../mesa.dat              
     rwf:   tmp/mesa.rwf                  rint:   tmp/mesa.rint            
     int:   tmp/mesa.int                  tint:   tmp/mesa.tint            
     gint:  tmp/mesa.gint                 rdint:  tmp/mesa.rdint           
     dint:  tmp/mesa.dint                 zint:   tmp/mesa.zint            
     ham:   tmp/mesa.ham                  moden:  tmp/mesa.moden           
     aoden: tmp/mesa.aoden                saoden: tmp/mesa.saoden          
     gden:  tmp/mesa.gden                 fci:    tmp/mesa.fci             

     machine size:        100000000
 kohn and bec files:
     kohn:  tmp/mesa.kohn                 kohndt: tmp/mesa.kohndt          
     grid:  tmp/mesa.grid                 orbs:   tmp/mesa.orbs            
     vstat: tmp/mesa.vstat                ylms:   tmp/mesa.ylms            
     bessel:tmp/mesa.bessel               knints: tmp/mesa.knints          
     tmat:  tmp/mesa.tmat                 blktmat:tmp/mesa.blktmt          
     optint:tmp/mesa.optint               atomci: tmp/mesa.atomci          
     lamdat:tmp/mesa.lamdat               bec:    tmp/mesa.bec             

     user defined maxsiz:  10000000
 title:
     DVR Code
 route:
     62//94;
     20//01;
 options:
     xprint=(m6294=(points,polynomials,matrix-elements)) generate-matrix-elements tes
     t-function=sine number-of-dimensions=1 xconvert=(left-endpoint=0.d0,right-endpoi
     nt=5.d0) space-dimension-1=x


                    grid and lobatto basis code

 first call to memory in link =     ptwt

 link            =     ptwt
 words requested =         61
 words available =   80000000
 words gotten    =         61

 link            =     poly
 words requested =        601
 words available =   79999939
 words gotten    =        601

 link            =  scratch
 words requested =        201
 words available =   79999338
 words gotten    =        201

 free memory link =  scratch
 words released   =        201

 link            =      mat
 words requested =       1001
 words available =   79999338
 words gotten    =       1001

 number-of-regions =    2
 left boundary condition = 0
 right boundary condition = 1

 edges =  0.00000000E+00  0.50000000E+01  0.15000000E+02

 link            =     hmat
 words requested =        649
 words available =   79998337
 words gotten    =        649

 link            =     norm
 words requested =         73
 words available =   79997688
 words gotten    =         73

 link            =        v
 words requested =         37
 words available =   79997615
 words gotten    =         37

 potential type = exponential                     
    using atomic units

     size of global basis set =   18

     region =   1 number of functions =   9
     starting function =   2 ending function =  10
     global starting function =   1 global ending function =   9
     bridge function = T

     region =   2 number of functions =   9
     starting function =   2 ending function =  10
     global starting function =  10 global ending function =  18
     bridge function = F

 link            =     diag
 words requested =        685
 words available =   79997578
 words gotten    =        685
eigenvalues                                                                     
 col          1
        -0.01079428
         0.01711543
         0.09538820
         0.22109480
         0.39244346
         0.60877159
         0.87159698
         1.15168630
         1.51870575
         2.15920027
         2.25663845
         3.26126426
         4.77499773
         5.94333202
        10.29934699
        16.72020846
        16.96192119
        37.85961843

 link            =     ramp
 words requested =         37
 words available =   79996893
 words gotten    =         37

       energy         k value       r-matrix       tan phase       phase     
  0.10000000E-04 0.44721360E-02 0.63108223E+01-0.38886257E-01-0.38866674E-01
  0.10000000E-03 0.14142136E-01 0.63476747E+01-0.12322035E+00-0.12260234E+00
  0.10000000E-01 0.14142136E+00 0.14354556E+02-0.15859836E+01-0.10082349E+01
  0.50000000E-01 0.31622777E+00-0.12452815E-01 0.28588033E+02 0.15358309E+01
  0.10000000E+00 0.44721360E+00-0.15874002E+02 0.34122942E+01 0.12857205E+01
  0.20000000E+00 0.63245553E+00 0.32247581E+01 0.17549814E+01 0.10528738E+01
  0.30000000E+00 0.77459667E+00-0.26473245E-01 0.13322843E+01 0.92691736E+00
  0.50000000E+00 0.10000000E+01 0.71913510E-01 0.98877331E+00 0.77975319E+00
  0.80000000E+00 0.12649111E+01 0.79051940E+00 0.77812021E+00 0.66125649E+00
  0.10000000E+01 0.14142136E+01-0.13410049E+00 0.67031509E+00 0.59052418E+00
  0.15000000E+01 0.17320508E+01 0.11485749E+01 0.26280837E+00 0.25699680E+00
  0.20000000E+01 0.20000000E+01 0.17029325E+01-0.47133650E+00-0.44045501E+00
  0.50000000E+01 0.31622777E+01 0.63580035E-01-0.11233164E+00-0.11186270E+00

 free memory link =     diag
 words released   =        685

 free memory link =     ramp
 words released   =         37

 free memory link =     hmat
 words released   =        649

 free memory link =     norm
 words released   =         73

 free memory link =        v
 words released   =         37

 free memory link =      mat
 words released   =       1001

 free memory link =     ptwt
 words released   =         61

 free memory link =     poly
 words released   =        601
 summary:
     link   entries               charges          core use      
                          cpu       sys      io         
        1         1       0.1       0.0       0.0         0
     6294         1       0.1       0.0       0.0  80000000
             total:       0.1       0.0       0.0  80000000
