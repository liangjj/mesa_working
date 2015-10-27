*deck @(#)fnlfun.f	1.1 9/8/91
c***begin prologue     fnlfun
c***date written       910805   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           m6020, link 6020, spline
c***author             schneider b.(lanl)
c***source             m6020
c***purpose            construct complex kohn free wave
c***                   
c***description        a single, complex wave whose imaginary part
c***                   is a regular ricatti-bessel function for all
c***                   values of rho and whose real part is a ricatti-
c***                   neumann function asymptotically, is constructed from
c***                   a bessel function and a solution of the first born
c***                   iterate to an exponential potential. 
c
c***references       
c
c***routines called    iosys, util and mdutil
c***end prologue       fnlfun
      subroutine fnlfun(final,scale,psil0,psil1,x,jnint,ynint,jnintm,
     1                  ynintm,n)
      implicit integer (a-z)
      real *8 wron, psil0, psil1, jnint, ynint, x
      real *8 jnintm, ynintm, scale
      complex *16 final, ai, f1, f2, c1, c2
      dimension psil0(n), psil1(n), final(n), x(n)
c----------------------------------------------------------------------c
c                 the solution needed is:                              c
c                    soln = (i*j  - q ) / x                            c
c                               l    l                                 c
c            where q behaves asymptotically as a neumann function      c
c            but is regular at the origin                              c 
c----------------------------------------------------------------------c
      ai=cmplx(0.d+00,1.d+00)
      wron=psil0(n)*psil1(n-1)-psil0(n-1)*psil1(n)
      f1=( -ynint + ai*jnint )
      f2=( -ynintm + ai*jnintm )
      c1=( f1*psil1(n-1) - f2*psil1(n) )/wron
      c2=( f2*psil0(n) - f1*psil0(n-1) )/wron
      do 10 i=1,n
         final(i) = ( c1*psil0(i) + c2*psil1(i) ) / x(i)
   10 continue
      scale=1.d+00/c2
      return
      end
