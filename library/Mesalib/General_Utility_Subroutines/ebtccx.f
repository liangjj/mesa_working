*deck ebtccx
c***begin prologue     ebtccx
c***date written       880423   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           complex matrix multiply, transpose
c***author             schneider, barry (lanl)
c***source             mylib 
c***                                         t
c***purpose            matrix multiply  a = b * c
c***description        vectorized matrix multiply
c***                   for complex a and c
c***                  
c***references         none
c
c***routines called    czero(mylib)
c***end prologue       ebtccx
      subroutine ebtccx(a,na,b,nb,c,nc,ni,nk,nj)
      implicit integer (a-z)
      real *8 b
      complex *16 a, c
      dimension a(na,nj), b(nb,ni), c(nc,nj)
      call czero(a,na*nj)
      do 10 k=1,nj
         do 20 j=1,nk
            do 30 i=1,ni
               a(i,k)=a(i,k)+b(j,i)*c(j,k)
   30       continue
   20    continue
   10 continue
      return
      end
