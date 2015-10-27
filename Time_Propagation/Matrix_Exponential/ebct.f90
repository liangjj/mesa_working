!deck @(#)ebct.f	5.1  11/6/94
      subroutine ebct(a,b,c,ni,nk,nj)
!***begin prologue     ebct
!***date written       850601  (yymmdd)
!***revision date      yymmdd  (yymmdd)
!***keywords           matrix, multiply,
!***author             saxe, paul (lanl)
!***source             @(#)ebct.f	5.1   11/6/94
!***purpose                                               t
!                      vectorized matrix operation:  a=b*c .
!***description
!                      call ebct(a,b,c,ni,nk,nj)
!                        a       output matrix, (ni,nj).
!                        b       input matrix, (ni,nk).
!                        c       input matrix, (nj,nk).
!
!***references
!***routines called    sgmm(clams)
!***end prologue       ebct
      implicit none
!
      integer                       :: ni, nj, nk
      real*8, dimension(ni,nj)      :: a
      real*8, dimension(ni,nk)      :: b
      real*8, dimension(nj,nk)      :: c
      real*8                        :: zero=0.d0, one=1.d0
!
!
      call sgemm('n','t',ni,nj,nk,one,b,ni,c,nj,zero,a,ni)
!
      return
      end
