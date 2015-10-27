*deck apbtcc
c***begin prologue     apbtcc
c***date written       880423   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           complex matrix multiply
c***author             schneider, barry (lanl)
c***source             mylib
c***                                             t
c***purpose            matrix multiply  a = a + b * c
c***description        vectorized matrix multiply
c***                   for complex a and c
c***                  
c***references         none
c
c***routines called    czero(mylib)
c***end prologue       apbtcc 
      subroutine apbtcc(a,b,c,ni,nk,nj)
      implicit integer (a-z)
      real *8 b
      complex *16 a,c
      dimension a(ni,nj), b(nk,ni), c(nk,nj)
      do 10 k=1,nj
         do 20 j=1,nk
            do 30 i=1,ni
               a(i,k)=a(i,k)+b(j,i)*c(j,k)
   30       continue
   20    continue
   10 continue
      return
      end
