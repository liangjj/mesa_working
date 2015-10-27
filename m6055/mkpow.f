*deck @(#)mkpow.f	1.1 9/8/91
c***begin prologue     mkpow
c***date written       xxxxxx   (yymmdd)
c***revision date      890427   (yymmdd)
c***keywords           kohn, integrals
c***author             rescigno, tom (llnl)
c***source             m6005
c***purpose            fill array ipow
c***references
c
c***routines called    none
c***end prologue       mkpow
      subroutine mkpow(ipow,lmax)
      implicit integer (a-z)
      dimension ipow(0:lmax)
      data zero, one/ 0, 1/
      ipow(0)=zero
      do 10 i=1,lmax
         ipow(i)=one
   10 continue
      return
      end
