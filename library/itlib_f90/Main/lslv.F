*deck lslv.f
c***begin prologue     lslv
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
c***end prologue       lslv
      subroutine lslv(a,b,ipvt,n,m,dim)
      implicit integer (a-z)
      real*8 a, b
      dimension a(dim,*), b(dim,m), ipvt(*) 
      common/io/inp, iout
      call sgetrf(n,n,a,dim,ipvt,info)
      if(info.ne.0) then
         call lnkerr('error from linear solve routine')
      endif   
      do 10 i=1,m
         call sgetrs('n',n,m,a,dim,ipvt,b,dim,info)
 10   continue   
      return
      end       



