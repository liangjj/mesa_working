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
      subroutine lsolve(a,b,ipvt,n,m,dim)
      implicit integer (a-z)
      real*8 a, b
      dimension a(dim,*), b(dim,m), ipvt(*) 
      common/io/inp, iout
      call sgefa(a,dim,n,ipvt,info)
      if(info.ne.0) then
         call lnkerr('error from linear solve routine')
      endif   
      do 10 i=1,m
         call sgesl(a,dim,n,ipvt,b(1,i),0)
 10   continue   
      return
      end       
