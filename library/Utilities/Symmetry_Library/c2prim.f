*deck @(#)c2prim.f	5.1  11/6/94
      subroutine c2prim(t,naxis,axis,told)
c     
      implicit integer (a-z)
c     
      real*8 t(3,3,naxis),told(3,3,naxis)
      character*(*) axis
c     
c     generate the additional c2 axes.  
c     the principal rotation axis is passed in axis.  the default is
c     the z -axis, in which case the additional ones are generated as
c     c2(y),c2(x). if the principal axis is different, generate additional ones
c     in the order of the cyclic permutation.
      if (axis.eq.' '.or.axis.eq.'z') then
c        c2(y), then c2(x)
         do 2 i=1,naxis
            do 1 j=1,3
               t(1,j,i)=-told(1,j,i)
               t(2,j,i)= told(2,j,i)
               t(3,j,i)=-told(3,j,i)
    1       continue
    2    continue
      else if (axis.eq.'y') then
c        c2(x), then c2(z)
         do 4 i=1,naxis
            do 3 j=1,3
               t(1,j,i)= told(1,j,i)
               t(2,j,i)=-told(2,j,i)
               t(3,j,i)=-told(3,j,i)
    3       continue
    4    continue
      else if (axis.eq.'x') then
c        c2(z), then c2(y)
         do 6 i=1,naxis
            do 5 j=1,3
               t(1,j,i)=-told(1,j,i)
               t(2,j,i)=-told(2,j,i)
               t(3,j,i)= told(3,j,i)
    5       continue
    6    continue
      else
         call lnkerr('c2prim: unrecognized axis')
      end if
c     
c     
      return
      end
