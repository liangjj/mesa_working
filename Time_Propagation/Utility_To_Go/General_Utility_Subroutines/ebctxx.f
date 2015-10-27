*deck @(#)ebctxx.f	1.3  8/8/91
      subroutine ebctxx(a,b,c,ni,nk,nj,na,nb,nc)
c***begin prologue     ebctxx
c***date written       850601  (yymmdd)
c***revision date      yymmdd  (yymmdd)
c***keywords           matrix, multiply,
c***author             saxe, paul (lanl)
c***source             @(#)ebctxx.f	1.3   8/8/91
c***purpose                                               t
c                      vectorized matrix operation:  a=b*c .
c***description
c                      call ebctxx(a,na,b,nb,c,nc,ni,nk,nj)
c                        a       output matrix, (na,nj).
c                        b       input matrix, (nb,nk).
c                        c       input matrix, (nc,nk).
c
c***references
c***routines called    sgmm(clams)
c***end prologue       ebctxx
      implicit integer(a-z)
c
      real*8 a(na,nj),b(nb,nk),c(nc,nk)
      real*8 zero,one
c
      parameter (zero=0.0d+0,one=1.0d+0)
c
c     call sgmm(ni,nk,nj,b,ni,c,nj,a,ni,2,1)
c     call rzero(a,ni*nj)
c     call mxmb(b,1,ni,c,nj,1,a,1,ni,ni,nk,nj)
c
c     call mxma(b,1,ni,c,nj,1,a,1,ni,ni,nk,nj)
      call sgemm('n','t',ni,nj,nk,one,b,nb,c,nc,zero,a,na)
c
      return
      end
