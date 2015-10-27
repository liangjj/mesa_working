*deck @(#)kzero.f	1.1  11/30/90
      subroutine kzero(x,nsmall,nbf)
      real*8 x(nbf,nbf)
c
      do 1 i=nsmall+1,nbf
         do 2 j=nsmall+1,nbf
            x(i,j)=0.d0
  2      continue
  1   continue
c
      return
      end
