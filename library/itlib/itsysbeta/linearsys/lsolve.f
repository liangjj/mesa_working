*deck lsolve.f
c***begin prologue     lsolve
c***date written       960723   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           linear system solve
c***author             schneider, barry (nsf)
c***source             
c***purpose            driver for direct linear system solve.
c***                   
c***references         
c
c***routines called    
c***end prologue       lsolve
      subroutine lsolve(a,b,work,ipvt,n,m,dim,lwork)
      implicit integer (a-z)
      real*8 a, b, work
      dimension a(dim,*), b(dim,m), ipvt(*) 
      common/io/inp, iout
      call dsysv('u',n,m,a,dim,ipvt,b,dim,work,lwork,info)
      if(info.ne.0) then
         stop
      endif   
      return
      end       
