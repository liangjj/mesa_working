*deck @(#)gshmdt.f	5.1  11/6/94
      subroutine gshmdt(a,is,mdim,nc,npt,rmax,rmin,rlim)
c***begin prologue     gshmdt.f
c***date written       yymmdd  
c***revision date      11/6/94      
c
c***keywords           
c***author             binkley, et al.(g82)
c***source             @(#)gshmdt.f	5.1   11/6/94
c***purpose            
c***description
c     form an orthonormal set of basis vectors from a set of given
c     vectors by schmidt orthogonalization.  vectors shorter than
c     rmin, or longer than rmax, or contributing less than rlim to
c     a new basis vector are discarded.  is is loaded with a list
c     of the vectors whish have been kept.  nc is the number of
c     basis functions, and npt is the number of input points, whish
c     are stored consecutively in a.  on exit, npt contains the
c     number of points kept.  mdim is the allocated dimension of a.
c
c***references
c
c***routines called
c
c***end prologue       gshmdt.f
      implicit none
c     --- input variables -----
      integer mdim,nc,npt
      real*8 rmax,rmin,rlim
c     --- input arrays (unmodified) ---
c     --- input arrays (scratch) ---
c     --- output arrays ---
      integer is(*)
      real*8 a(mdim,*)
c     --- output variables ---
c     --- scratch arrays ---
c     --- local variables ---
      integer i,j,ii,k,i1
      real*8 r,rr
      real*8 sdot
      real*8 zero,one
c
      parameter (zero=0.0d+00,one=1.0d+00)
c
 1030 format(11(10x,10f12.6/))
 1031 format(40i3)
c
c
      i=1
      ii=0
   20 ii=ii+1
      is(i)=ii
      rr = sqrt(sdot(nc,a(1,i),1,a(1,i),1))
      if((npt.gt.1).and.(rr.lt.rmin)) goto 115
      if((i.gt.1).and.(rr.gt.rmax)) goto 115
      call smul(a(1,i),a(1,i),(one/rr),nc)
      rr=one
      i1=i-1
      if(i1.lt.1) goto 110
          do 70 k=1,i1
              r = sdot(nc,a(1,i),1,a(1,k),1)
   70         call vwxs(a(1,i),a(1,i),a(1,k),-r,+1,nc)
          rr = sqrt(sdot(nc,a(1,i),1,a(1,i),1))
          if(rr.lt.rlim) goto 115
          call smul(a(1,i),a(1,i),(one/rr),nc)
  110 if(i.ge.npt) goto 120
      i=i+1
      goto 20
c
c     --- discard any previous points that are unsuitable.
  115 npt=npt-1
      if(i.gt.npt) goto 120
      do 100 j=1,nc
          do 100 k=i,npt
  100         a(j,k)=a(j,k+1)
      goto 20
  120 continue
c
c 
      return
      end
