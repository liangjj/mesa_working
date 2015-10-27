*deck matcpy.f
c***begin prologue     matcpy
c***date written                (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           dimensioned matrix copy
c***author             schneider, barry (nsf)
c***source             
c***purpose            
c***                   
c***                                                          
c***references         
c
c***routines called    
c***end prologue       matcpy
      subroutine matcpy(a,b,n,m,ni,nj)
      implicit integer (a-z)
      real*8 a, b
      dimension a(ni,*), b(nj,*)
      common/io/inp, iout
      do 10 i=1,n
         call copy(a(1,i),b(1,i),m)
 10   continue
      return   
      end       
