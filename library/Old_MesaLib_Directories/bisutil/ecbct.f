*deck ecbct
c***begin prologue     ecbct
c***date written       880423   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           complex matrix multiply
c***author             schneider, barry (lanl)
c***source             mylib
c***                                             t
c***purpose            matrix multiply  a = b * c
c***description        vectorized matrix multiply 
c***                   for complex a and b
c***                  
c***references         none
c
c***routines called    czero(mylib)
c***end prologue       ecbct
      subroutine ecbct(a,b,c,ni,nk,nj)
      implicit integer (a-z)
      real *8 c
      complex *16 a, b
      dimension a(ni,nj), b(ni,nk), c(nj,nk)
      call czero(a,ni*nj)
      do 10 k=1,nj
         do 20 j=1,nk
            do 30 i=1,ni
               a(i,k)=a(i,k)+b(i,j)*c(k,j)
   30       continue
   20    continue
   10 continue
      return
      end
