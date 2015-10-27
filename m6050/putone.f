*deck @(#)putone.f
c***begin prologue     putone
c***date written       920501   (yymmdd)
c***revision date               (yymmdd)
c***keywords           
c***author             schneider, barry (nsf)
c***source             m6050
c***purpose            unit matrix
c***references         none
c
c***routines called    none
c***end prologue       grdprn
      subroutine putone(e,ev,n,m)
      implicit integer (a-z)
      complex *16 e, ev, one
      dimension e(m), ev(n,m)
      data one /(1.d0,1.d0)/
      call czero(e,m)
      call czero(ev,n*m)
      do 10 i=1,m
         e(i)=one
         do 20 j=1,n
            ev(j,i)=one
   20    continue 
   10 continue
      return
      end
