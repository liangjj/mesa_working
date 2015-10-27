*deck chk.f
c***begin prologue     chk
c***date written       951229   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           chebychev polynomials
c***author             schneider, barry (nsf)
c***source             math
c***purpose            check orthonormality of polynomials.
c***                   
c***references         
c
c***routines called    
c***end prologue       chk
      subroutine chk(p,wts,scr,dum,n,npts)
      implicit integer (a-z)
      real*8 p, wts, scr, dum
      character*80 title
      dimension p(npts,0:n-1), wts(npts), scr(*), dum(*)
      common/io/inp, iout
      call vsqrt(dum,wts,npts)
      call copy(p,scr,npts*n)
      call vmmul(dum,scr,scr,npts,n) 
      call ebtc(dum,scr,scr,n,npts,n)
      title='overlap matrix'
      call prntrm(title,dum,n,n,n,n,iout)
      call ebct(dum,scr,scr,npts,n,npts)
      title='completeness matrix'
      call prntrm(title,dum,npts,npts,npts,npts,iout)
      return
      end       
