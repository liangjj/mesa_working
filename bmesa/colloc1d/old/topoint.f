*deck topoint.f
c***begin prologue     topoint
c***date written       950704   (yymmdd)
c***revision date               (yymmdd)
c***keywords           functions, inverse
c***author             schneider, barry(nsf)
c***source             
c***purpose            calculate inverse needed to convert from spectral
c***                   to grid representation.
c***
c***
c***references
c
c***routines called
c***end prologue
      subroutine topoint(a,work,ipvt,n,prnt)
      implicit integer (a-z)
      dimension a(n,n), work(n), ipvt(n), det(2)
      real*8 a, work, det
      character*80 title
      logical prnt
      common /io/ inp, iout
      call sgefa(a,n,n,ipvt,info)
      if (info.ne.0) then
          call lnkerr('error in inverse collocation matrix')
      else          
          call sgedi(a,n,n,ipvt,det,work,1)
      endif
      if (prnt) then
          title='inverse of collocation matrix'
          call prntrm(title,a,n,n,n,n,iout)          
      endif
      return
      end
