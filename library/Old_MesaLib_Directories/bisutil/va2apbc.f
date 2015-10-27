*deck va2apbc
c***begin prologue     va2apbc
c***date written       980923   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           vector add plus multiply
c***author             schneider, barry (lanl)
c***source             utila
c***purpose            vector add plus multiply  va = va + vb * vc
c***description        vectorized vector multiply
c***references         none
c
c***routines called    
c***end prologue       va2apbc 
      subroutine va2apbc(va,vb,vc,n,m)
      implicit integer (a-z)
      real*8 va, vb, vc
      dimension va(n,m), vb(n,m), vc(n,m)
      do 10 i=1,m
         do 20 j=1,n
            va(j,i) = va(j,i) + vb(j,i)*vc(j,i)
 20      continue   
 10   continue
      return
      end
