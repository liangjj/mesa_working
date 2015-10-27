*deck ecbtc
c***begin prologue     ecbtc
c***date written       880423   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           complex matrix multiply, transpose
c***author             schneider, barry (lanl)
c***source             mylib 
c***                                         t
c***purpose            matrix multiply  a = b * c
c***description        vectorized matrix multiply
c***                   for complex a and b
c***                  
c***references         none
c
c***routines called    czero(mylib)
c***end prologue       ecbtc
      subroutine ecbtc(a,b,c,ni,nk,nj)
      implicit integer (a-z)
      real *8 c
      complex *16 a, b
      dimension a(ni,nj), b(nk,ni), c(nk,nj)
      call czero(a,ni*nj)
      do 10 k=1,nj
         do 20 j=1,nk
            do 30 i=1,ni
               a(i,k)=a(i,k)+b(j,i)*c(j,k)
   30       continue
   20    continue
   10 continue
      return
      end
