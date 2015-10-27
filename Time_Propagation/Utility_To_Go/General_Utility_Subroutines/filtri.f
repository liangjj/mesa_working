*deck filtri
c***begin prologue     filtri
c***date written       920209   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           symmetrize matrix
c***author             schneider, barry (lanl)
c***source             mylib 
c***                                         
c***purpose            move lower triangle of symmetric matrix
c***                   into upper triangle
c***                   
c***                  
c***references         none
c
c***routines called    
c***end prologue       filtri
      subroutine filtri(a,n)
      implicit integer (a-z)
      real *8 a
      dimension a(n,n)
      do 10 i=1,n
         do 20 j=1,i
            a(j,i)=a(i,j)
   20    continue
   10 continue
      return
      end
