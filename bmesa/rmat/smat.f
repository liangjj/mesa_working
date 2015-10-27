*deck smat
c***begin prologue     smat
c***date written       940213   (yymmdd)
c***revision date               (yymmdd)
c***keywords           bessel, overlap
c***author             schneider, barry (nsf)
c***source             
c***purpose            overlap of numerical basis functions
c***description        
c***references       
c
c***routines called
c***end prologue       smat
      subroutine smat(f,s,npts,n,prnt)
      implicit integer (a-z)
      common /io/ inp, iout
      real*8 f, s
      character*80 title
      logical prnt
      dimension f(npts,n), s(n,n)
c     calculate overlap integral
      call ebtc(s,f,f,n,npts,n)
      if (prnt) then
          title='primitive overlap matrix'
          call prntrm(title,s,n,n,n,n,iout)
      endif
      return
      end
