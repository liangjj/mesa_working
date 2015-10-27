*deck @(#)printr.f	1.1  11/30/90
      subroutine printr(x,leni,lenj)
c
c***begin prologue
c***date written       871022   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c
c***keywords
c***author             lengsfield, byron (brl)
c***source             @(#)printr.f	1.1   11/30/90
c
c***purpose
c
c***description
c
c***references
c
c***routines called    (none)
c
c***end prologue
c
      implicit real*8(a-h,o-z)
c
      common /io/ inp,iout
c
cc
cmp   extended dummy x
cc
      dimension x(leni,lenj)
      write(iout,100)
      je=0
      do 11 k=1,lenj,4
         js=je+1
         je=je+4
         je=min(je,lenj)
         write(iout,13)(ii,ii=js,je)
         do 10 i=1,leni
            write(iout,12) i,(x(i,j),j=js,je)
 10      continue
 11   continue
 12   format(2x,i4,4(2x,f14.7))
 13   format(6x,4(7x,i2,7x))
c
 100  format(//,'  rectangular matrix  ',/)
c
      return
      end
