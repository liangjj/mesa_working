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
      subroutine lsolve(a,ca,b,cb,ipvt,n,m,dim,mattyp)
      implicit integer (a-z)
      character*(*) mattyp
      real*8 a, b
      complex*16 ca, cb
      dimension a(dim,*), ca(dim,*), b(dim,m), cb(dim,m), ipvt(*) 
      common/io/inp, iout
      if(mattyp.eq.'complex') then 
        call cgefa(ca,dim,n,ipvt,info)
        if(info.ne.0) then
           call lnkerr('error from linear solve routine')
        endif
        do 10 i=1,m
           call cgesl(ca,dim,n,ipvt,cb(1,i),0)
 10     continue   
      else
        call sgefa(a,dim,n,ipvt,info)
        if(info.ne.0) then
           call lnkerr('error from linear solve routine')
        endif   
        do 20 i=1,m
           call sgesl(a,dim,n,ipvt,b(1,i),0)
 20     continue   
      endif
      return
      end       
