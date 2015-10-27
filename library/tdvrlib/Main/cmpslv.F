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
      subroutine cmpslv(ham,rhs,ipvt,n,m,prn)
      implicit integer (a-z)
      real*8 ham, rhs
      character*80 title
      logical prn
      dimension  ham(n,n), rhs(n,m), ipvt(n) 
      common/io/inp, iout
      write(iout,1) n, m
      call sgetrf(n,n,ham,n,ipvt,info)
      call sgetrs('n',n,m,ham,n,ipvt,rhs,n,info)
      if(prn) then
         title='time-dependent coefficients'
         call prntrm(title,rhs,n,m,n,m,iout)
      endif
      return
 1    format(/,5x,'direct solution of complex linear system of:',/,5x,
     1 '                                            size = ',i4,/,5x,
     2 '                                            no.rhs = ',i3) 
      end       
