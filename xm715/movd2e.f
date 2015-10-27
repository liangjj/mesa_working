*deck %W%  %G%
      subroutine movd2e(ld2e,d2e,nd2e,dercen,ndcen,npass)
c***begin prologue     movd2e.f
c***date written       840427   (yymmdd)
c***revision date      11/6/94
c
c***keywords
c***author             saxe, paul (lanl)
c***source             %W%   %G%
c
c***purpose            to invoke translational invariance and move a
c   shell-block contribution of second derivatives to the total matrix.
c
c***description
c
c***references
c
c***routines called    (none)
c
c***end prologue       movd2e.f
      implicit none
c     --- input variables -----
      integer nd2e,ndcen,npass
c     --- input arrays (unmodified) ---
      integer dercen(4)
      real*8 ld2e(78)
      real*8 d2e(nd2e)
c     --- input arrays (scratch) ---
c     --- output arrays ---
c     --- output variables ---
c     --- scratch arrays ---
c     --- local variables ---
      integer atom(4)
      integer inp,iout
      integer iadd
      integer i,j,l,last,coord1,coord2
      integer ll,jj,lj,c1,c2,ic,jc,ii,ij,n2
      integer top
c
      common /io/ inp,iout
c
      iadd(i,j)=i*(i-1)/2+j
c
c     --- use translational invariance to fill out the list ---
c           xlxi=-xixi-xjxi-...
      last=ndcen+1
      l=3*(last-1)
      do 4 coord1=1,3
         ll=iadd(l+coord1,0)
         do 3 c2=1,ndcen
            j=3*(c2-1)
            do 2 coord2=1,3
               jc=j+coord2
               jj=iadd(jc,0)
               lj=ll+jc
               ld2e(lj)=0.0d+00
               do 1 c1=1,ndcen
                  ic=3*(c1-1)+coord1
                  if (ic.le.jc) then
                     ld2e(lj)=ld2e(lj)-ld2e(jj+ic)
                  else
                     ld2e(lj)=ld2e(lj)-ld2e(iadd(ic,jc))
                  end if
 1             continue
 2          continue
 3       continue
 4    continue
c
c     --- and xlxl=-xlxi-xlxj-.... ---
      do 7 coord1=1,3
         ll=iadd(l+coord1,0)
         do 6 coord2=1,coord1
            lj=ll+l+coord2
            ld2e(lj)=0.0d+00
            do 5 c2=1,ndcen
               ic=3*(c2-1)+coord2
               ld2e(lj)=ld2e(lj)-ld2e(ll+ic)
 5          continue
 6       continue
 7    continue
c
c     --- set up the order of atoms ---
      atom(1)=dercen(1)
      if (npass.eq.1) then
         atom(2)=dercen(2)
         atom(3)=dercen(3)
         atom(4)=dercen(4)
      else if (npass.eq.2) then
         atom(2)=dercen(2)
         atom(3)=dercen(3)
      else if (npass.eq.3) then
         atom(2)=dercen(3)
         atom(3)=dercen(2)
      else if (npass.eq.4) then
         atom(2)=dercen(2)
      end if
c
c     --- need to double eg. xixj if i=j ---
      if (atom(1).eq.atom(2)) then
         ld2e(7)=ld2e(7)*2.0d+00
         ld2e(12)=ld2e(12)*2.0d+00
         ld2e(18)=ld2e(18)*2.0d+00
      end if
c
c     --- transfer to the final, full second-derivative matrix
      n2=0
      do 11 c1=1,last
         i=3*(atom(c1)-1)
         do 10 coord1=1,3
            ic=i+coord1
            ii=iadd(ic,0)
            do 9 c2=1,c1
               j=3*(atom(c2)-1)
               if (c2.eq.c1) then
                  top=coord1
               else
                  top=3
               end if
               do 8 coord2=1,top
                  jc=j+coord2
                  if (ic.ge.jc) then
                     ij=ii+jc
                  else
                     ij=iadd(jc,ic)
                  end if
                  n2=n2+1
                  d2e(ij)=d2e(ij)+ld2e(n2)
 8             continue
 9          continue
 10      continue
 11   continue
c
c
      return
      end
