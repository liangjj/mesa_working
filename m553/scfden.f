*deck @(#)scfden.f	5.1  11/6/94
      subroutine scfden(cr,nobt,n12)
c
c***begin prologue     scfden
c***date written       871022   (yymmdd)
c***revision date      871118   (yymmdd)
c   18 november 1987   pws and bhl lanl
c   2-particle dm  dm(ii,jj) corrected
c
c***keywords
c***author             lengsfield, byron (brl)
c***source             @(#)scfden.f	5.1   11/6/94
c
c***purpose
c
c***description
c
c***references
c
c***routines called    (none)
c
c***end prologue       scfden
c
      implicit real*8(a-h,o-z)
      dimension cr(n12,n12)
c
      common /io/ inp,iout
c
      ioff(i,j)=i*(i-1)/2+j
c
c
      call rzero(cr,n12)
      do 1 i=1,nobt
         cr(ioff(i,i),1)=1.0d+00
 1    continue
c
c      ix=0
c      do 21 i=1,nobt
c         ix=ix+i
c         cr(ix)=1.0d+00
c 21   continue
c
      call iosys('write real "mo 1pdm" on mcscr',n12,cr,0,' ')
c
      call rzero(cr,n12**2)
c
      do 3 i=1,nobt
         ii=ioff(i,i)
         do 2 j=1,nobt
            if (j.ne.i) then
               jj=ioff(j,j)
               cr(ii,jj)=0.5d+00
            end if
 2       continue
 3    continue
c
      do 5 i=2,nobt
         ii=ioff(i,0)
         do 4 j=1,i-1
            ij=ii+j
            cr(ij,ij)=-0.25d+00
 4       continue
 5    continue
c
c$$$      call rzero(cr,ntot)
c$$$      ii=0
c$$$      do 22 i=1,nobt
c$$$         ii=ii+i
c$$$         ix=(ii*(ii-1))/2
c$$$         jj=0
c$$$         do 23 j=1,i
c$$$            jj=jj+j
c$$$            cr(ix+jj)=1.d0
c$$$ 23      continue
c$$$         cr(ix+jj)=0.d0
c$$$ 22   continue
c$$$c
c$$$      if(nobt.gt.1) then
c$$$c
c$$$c    open-shell exchange density
c$$$c
c$$$         ii=1
c$$$         do 24 i=2,nobt
c$$$            i1=i-1
c$$$            do 25 j=1,i1
c$$$               ij=ii+j
c$$$               index=(ij*(ij+1))/2
c$$$               cr(index)=-1.0d+00
c$$$ 25         continue
c$$$            ii=ii+i
c$$$ 24      continue
c$$$c
c$$$      endif
c
      call iosys('write real "mo 2pdm" on mcscr',n12**2,cr,0,' ')
c
      return
      end
