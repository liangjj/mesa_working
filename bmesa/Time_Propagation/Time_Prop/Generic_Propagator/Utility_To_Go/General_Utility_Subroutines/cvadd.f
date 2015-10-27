*deck cvadd
c***begin prologue     cvadd
c***date written       890602   (yymmdd)
c***revision date               (yymmdd)
c***keywords           complex matrix add
c***author             schneider, barry (lanl)
c***source             mylib
c***purpose            add two complex matrices
c***          
c***description        
c***                   
c***                   
c
c***references       
c
c***routines called   
c***end prologue      
      subroutine cvadd(a,n1,b,n2,n3,n4)
      implicit integer (a-z)
      complex*16 a, b
      dimension a(n1,n4), b(n2,n4)
      do 10 i=1,n3
         do 20 j=1,n4
            a(i,j)=a(i,j)+b(i,j)
   20    continue
   10 continue
      return
      end
