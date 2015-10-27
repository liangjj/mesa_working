*deck @(#)scdm1.f	1.1  11/30/90
      subroutine scdm1(dm,nao,nij)
      implicit integer(a-z)
      real*8 dm(nij)
c
      ix=0
      do 10 i=1,nao
         ix=ix+i
         dm(ix)=dm(ix)*.5d+00
  10  continue
c
      do 20 i=1,nij
         dm(i)=dm(i)*2.d+00
   20 continue
c
      return
      end
