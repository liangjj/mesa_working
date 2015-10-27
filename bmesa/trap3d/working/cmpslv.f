*deck cmpslv.f
c***begin prologue     cmpslv
c***date written       960723   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           complex linear equations
c***author             schneider, barry (nsf)
c***source             
c***purpose            driver for complex linear system solve.
c***                   
c***references         
c
c***routines called    
c***end prologue       cmpslv
      subroutine cmpslv(ham,rhs,ipvt,n,m)
      implicit integer (a-z)
      complex*16 ham, rhs
      character*80 title
      dimension  ham(n,n), rhs(n,m), ipvt(n) 
      common/io/inp, iout
      call cgefa(ham,n,n,ipvt,info)
      call cgesl(ham,n,n,ipvt,rhs,0)
      title='time-dependent coefficients'
      call prntcm(title,rhs,n,m,n,m,iout)
      return
      end       
