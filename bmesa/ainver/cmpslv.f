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
      subroutine cmpslv(ham,rhs,ipvt,n,m,prnt)
      implicit integer (a-z)
      real*8 ham, rhs
      character*80 title
      logical prnt
      dimension  ham(n,n), rhs(n,m), ipvt(n) 
      common/io/inp, iout
      write(iout,1) n, m
      call sgefa(ham,n,n,ipvt,info)
      call sgesl(ham,n,n,ipvt,rhs,0)
      if(prnt) then
         title='time-dependent coefficients'
         call prntrm(title,rhs,n,m,n,m,iout)
      endif
      return
 1    format(/,5x,'direct solution of complex linear system of:',/,5x,
     1 '                                            size = ',i4,/,5x,
     2 '                                            no.rhs = ',i3) 
      end       
