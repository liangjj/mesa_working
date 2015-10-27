*deck cvconj
c***begin prologue     cvconj
c***date written       890602   (yymmdd)
c***revision date               (yymmdd)
c***keywords           complex conjugate
c***author             schneider, barry (lanl)
c***source             mylib
c***purpose           
c***          
c***description        
c***                   
c***                   
c
c***references       
c
c***routines called   
c***end prologue      
      subroutine cvconj(a,b,n)
      implicit integer (a-z)
      complex*16 a, b
      dimension a(n), b(n)
      do 10 i=1,n
         b(i) = conjg(a(i))
   10 continue
      return
      end
