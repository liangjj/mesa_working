*deck @(#)rzero.f	5.1  11/6/94
      subroutine rzero(a,n)
c***begin prologue     rzero
c***date written       850601  (yymmdd)
c***revision date      yymmdd  (yymmdd)
c***keywords           vector, initialize
c***author             saxe, paul (lanl)
c***source             @(#)rzero.f	5.1   11/6/94
c***purpose            vectorized vector initialization: a=0.0 .
c***description
c                      call rzero(a,n)
c                        a         output vector, declared real.
c                        n         vector length.
c
c***references
c***routines called    (none)
c***end prologue       rzero
      real*8 a(n)
c
      do 1 i=1,n
           a(i)=0.0d+00
    1 continue
c
      return
      end
