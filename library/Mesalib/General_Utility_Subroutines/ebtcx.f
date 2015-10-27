*deck ebtcx
c***begin prologue     ebtcx
c***date written       880423   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           matrix multiply, transpose
c***author             schneider, barry (lanl)
c***source             mylib 
c***                                         t
c***purpose            matrix multiply  a = b * c
c***description        vectorized matrix multiply
c***                   
c***                  
c***references         none
c
c***routines called    sgemm
c***end prologue       ebtcx
      subroutine ebtcx(a,na,b,nb,c,nc,ni,nk,nj)
      implicit integer (a-z)
      real *8 a, b, c
      real *8 zero, one
      dimension a(na,nj), b(nb,ni), c(nc,nj)
      parameter (zero=0.0d+0,one=1.0d+0)
      call sgemm('t','n',ni,nj,nk,one,b,nb,c,nc,zero,a,na)
c      call sgmm(ni,nk,nj,b,nb,c,nc,a,na,4,1)
      return
      end
