
*
*     ------------------------------------------------------------------
*               i n i t
*     ------------------------------------------------------------------
*
*       initializes basic constants of the program including those
*   which define the average energy of a configuration.
*
*
      SUBROUTINE init
      IMPLICIT DOUBLE PRECISION(a-h,o-z)
      PARAMETER (nt=100,kx=8,nwd=10)
      COMMON /eav/cca(10),ccb(35)
      COMMON /fact/gam(100)
*
      COMMON /PARAM/ z,h,t(nt),fine,mass,ks,ns,nwf,nclosd,ncfg,ib,ic,
     :               id,nscf,no,kmax,lmax,noff
*
*  ***** set the factorial needed by rme
*
      call factrl(32)
*
*  *****  set fine structure constant
*
      fine = 0.25d0/(137.036)**2
*
*  *****  set the starting point, step size, and related parameters
*
      rho = -4.d0
      h   = 1./16.d0
      h1 = h/1.5d0
      h3 = h/3.d0
      ch = h*h/12.d0
      eh = dexp(-h)
*
*  *****  average interactions for equivalent electrons
*
*  *****  p - p
*
      cca(1) = 2.d0/25.d0
*
*  *****  d - d
*
      cca(2) = 2.d0/63.d0
      cca(3) = 2.d0/63.d0
*
*  *****  f - f
*
      cca(4) =   4.d0/ 195.d0
      cca(5) =   2.d0/ 143.d0
      cca(6) = 100.d0/5577.d0
*
*  *****  g - g
*
      cca(7) =   20.d0/  1309.d0
      cca(8) =  162.d0/ 17017.d0
      cca(9) =   20.d0/  2431.d0
      cca(10) = 4410.d0/371943.d0
*
*
*  ***** average interactions for non-equivalent electrons
*
*  *****  s - ( s, p, d, f, g )
*
      ccb(1) = 1.d0/ 2.d0
      ccb(2) = 1.d0/ 6.d0
      ccb(3) = 1.d0/10.d0
      ccb(4) = 1.d0/14.d0
      ccb(5) = 1.d0/18.d0
*
*  *****  p - ( p, d, f, g )
*
      ccb(6) = 1.d0/  6.d0
      ccb(7) = 1.d0/ 15.d0
      ccb(8) = 1.d0/ 15.d0
      ccb(9) = 3.d0/ 70.d0
      ccb(10) = 3.d0/ 70.d0
      ccb(11) = 2.d0/ 63.d0
      ccb(12) = 2.d0/ 63.d0
      ccb(13) = 5.d0/198.d0
*
*  *****  d - ( d, f, g )
*
      ccb(14) =  1.d0/ 10.d0
      ccb(15) =  1.d0/ 35.d0
      ccb(16) =  1.d0/ 35.d0
      ccb(17) =  3.d0/ 70.d0
      ccb(18) =  2.d0/105.d0
      ccb(19) =  5.d0/231.d0
      ccb(20) =  1.d0/ 35.d0
      ccb(21) = 10.d0/693.d0
      ccb(22) =  5.d0/286.d0
*
*  *****  f - ( f, g )
*
      ccb(23) =  1.d0/  14.d0
      ccb(24) =  2.d0/ 105.d0
      ccb(25) =  1.d0/  77.d0
      ccb(26) = 50.d0/3003.d0
      ccb(27) =  2.d0/  63.d0
      ccb(28) =  1.d0/  77.d0
      ccb(29) = 10.d0/1001.d0
      ccb(20) = 35.d0/2574.d0
*
*  *****  g - ( g )
*
      ccb(31) =   1.d0/   18.d0
      ccb(32) =  10.d0/  693.d0
      ccb(33) =   9.d0/ 1001.d0
      ccb(34) =  10.d0/ 1287.d0
      ccb(35) = 245.d0/21879.d0
      return
      end
