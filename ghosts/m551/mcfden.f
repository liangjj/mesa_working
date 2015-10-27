*deck %W%  %G%
      subroutine mcfden(den,temp,nao)
c
c***begin prologue
c***date written       871022   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c
c***keywords
c***author             lengsfield, byron (brl)
c***source             %W%   %G%
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
c      extended dummy den,temp
c
      dimension den(2),temp(2)
c
      ntot=(nao*(nao+1))/2
c
      do 10 i=1,ntot
         temp(i)=-den(i)
 10   continue
c
      ix=0
      do 20 i=1,nao
         ix=ix+i
         temp(ix)=2.d0*temp(ix)
 20   continue
c
      return
      end
