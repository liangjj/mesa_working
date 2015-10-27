*deck @(#)minvrt.f	5.1  11/6/94
      subroutine minvrt(a,lda,n,det,ipvt,work)
c***begin prologue     minvrt.f
c***date written       850601  (yymmdd)
c***revision date      11/6/94
c   august 2, 1991     rlm at lanl
c      generalizing to handle other cases.
c***keywords           matrix, invert
c***author             martin, richard (lanl)
c***source             @(#)minvrt.f	5.1   11/6/94
c***purpose            vectorized matrix inversion.
c***description
c
c     call minvrt(a,lda,n,det,ipvt,work)
c
c     invert the real symmetric matrix a.
c     a      ... the real symmetric input matrix, real(lda,n).
c     lda    ... the leading dimension of the matrix a.
c     n      ... the order of the matrix a.
c     det    ... the determinant of a.
c     ipvt   ... scratch space, integer(n).
c     work   ... scratch space, real(n).
c
c***references
c***routines called    sgefa(clams), sgedi(clams), lnkerr(mdutil)
c***end prologue       minvrt.f
      implicit integer(a-z)
      real*8 a(lda,n),dettmp(2),work(n),det,rcond
      real*8 zero
      integer ipvt(n)
c
      parameter (zero=0.0d+0)
c
c     factor the matrix by gaussian elimination.
      call sgeco(a,lda,n,ipvt,rcond,work)
      if (rcond.eq.zero) then
         call lnkerr('matrix singular in minvrt')
      endif
c
c
      call sgedi(a,lda,n,ipvt,dettmp,work,11)
      det=dettmp(1)*10**dettmp(2)
c
c
      return
      end
