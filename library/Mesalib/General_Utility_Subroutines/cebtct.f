*deck @(#)cebtct.f	1.3  8/8/91
      subroutine cebtct(a,b,c,ni,nk,nj)
c***begin prologue     cebtct
c***date written       910903  (yymmdd)
c***revision date      yymmdd  (yymmdd)
c***keywords           matrix, multiply,complex
c***author             schneider, barry (lanl)
c***source             @(#)cebtct.f	1.3   8/8/91
c***purpose                                             t  t
c                      vectorized matrix operation:  a=b *c .
c***description
c                      call ebtct(a,b,c,ni,nk,nj)
c                        a       complex output matrix, (ni,nj).
c                        b       complex input matrix, (nk,ni).
c                        c       complex input matrix, (nj,nk).
c
c***references
c***routines called    sgmm(clams)
c***end prologue       cebtct
      implicit  integer(a-z)
c
      complex *16 a(ni,nj),b(nk,ni),c(nj,nk)
      complex *16 zero, one
c
      zero=cmplx(0.d0,0.d0)
      one=cmplx(1.d0,0.d0)
c
      call cgemm('t','t',ni,nj,nk,one,b,nk,c,nj,zero,a,ni)
c
      return
      end
