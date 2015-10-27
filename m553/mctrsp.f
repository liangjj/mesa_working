*deck @(#)mctrsp.f	5.1  11/6/94
      subroutine mctrsp(xk,xkt,nbr,nbs)
c
c***begin prologue
c***date written       871022   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c
c***keywords
c***author             lengsfield, byron (brl)
c***source             @(#)mctrsp.f	5.1   11/6/94
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
cc
cmp   extended dummy xk,xkt
cc
      dimension xk(nbr,nbs),xkt(nbs,nbr)
      common / number / zero,pt5,one,two,four,eight
c
c------------------------------------------------------
c  this program is used to transpose exchange operators
c------------------------------------------------------
c
      do 20 i=1,nbs
         do 10 j=1,nbr
            xkt(i,j)=xk(j,i)
 10      continue
 20   continue
c
      return
      end
