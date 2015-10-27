*deck @(#)clvcmp.f	1.1 9/8/91
c***begin prologue     clvcmp
c***date written       880423   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           clvcmp, link 1106, kohn variational
c***author             schneider, barry (lanl), rescigno, tom(llnl)  
c***source             m1106
c***purpose            solve complex linear systems
c***description        solve a set(s) of linear equations in which
c***                   the matrix and the right hand side are
c***                   complex.
c***                  
c***references         schneider and rescigno, physical review
c
c***routines called    iosys, util and mdutil
c***end prologue       clvcmp
      subroutine clvcmp (hamcc,ipvt,rhsf,n,m,ifac,prnt,fil)
      implicit integer(a-z)
      complex*16  hamcc, rhsf
      character *(*) ifac
      character *24 fil
      character *80 title
      logical prnt
      dimension hamcc(n,n), rhsf(n,m)
      dimension ipvt(n)
      common /io/ inp, iout
c----------------------------------------------------------------------c
c     compute lu factorization of complex matrix                       c
c----------------------------------------------------------------------c
      if (ifac.eq.'factor') then
         if (fil.ne.'no write') then
             call iosys ('write real '//fil//' to tmat',2*n*n,
     1                    hamcc,0,' ')
         endif
         call cgefa (hamcc,n,n,ipvt,info)
      else
c----------------------------------------------------------------------c
c                solve equations                                       c
c----------------------------------------------------------------------c
         if (fil.ne.'no write') then
             call iosys ('write real '//fil//' to tmat',2*m*n,
     1                    rhsf,0,' ')
         endif
         do 20 i=1,m
            call cgesl (hamcc,n,n,ipvt,rhsf(1,i),0)
   20    continue
         if (prnt) then
             title='complex-solution'
             call prntcm(title,rhsf,n,m,n,m,iout)
         endif
      endif
      return
      end
