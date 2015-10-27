*deck %W%  %G%
      subroutine zerodm(den,nao,iobcas)
c
c***begin prologue
c***date written       880128   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c
c***keywords
c***author             lengsfield, byron (brl)
c***source             %W%   %G%
c
c***purpose       zero 1-particle density matrix elements
c                 between invariant orbital sets so the
c                 outgoing molecular orbitals remain
c                 variational..orbitals within an invariant set
c                 are rotated so they diagonalize a sub-block
c                 of the 1-particle density matrix
c
c***description
c
c***references
c
c***routines called    (none)
c
c***end prologue
      implicit real*8(a-h,o-z)
      logical debug
      real*8 den(*)
      dimension iobcas(*)
c
      data debug/.false./
      common /io/ inp,iout
c
      iflag=-1
c
      if(nao.le.1) return
ccc
c ---- test the number of orbital sets
ccc
      iset=iobcas(1)
      do 1 i=2,nao
         if(iobcas(i).eq.iset) go to 2
         iflag=0
   2  continue
   1  continue
ccc
c ---- if there is only one set return
ccc
      if(iflag.lt.0) return
c
      if (debug) then
c        naot=nao*(nao+1)/2
c        write(iout,*)' incoming density matrix '
c        call print(den,naot,nao,iout)
      end if
c
ccc
c ---- zero matrix elements coupling invariant sets
ccc
      ix=0
      do 3 i=1,nao
         iset=iobcas(i)
         do 4 j=1,i
            jset=iobcas(j)
            ix=ix+1
            if(iset.ne.jset) then
               iflag=1
               den(ix)=0.d0
            end if
  4      continue
  3   continue
c
      if(iflag.ne.0) then
         write(iout,5)
  5      format(/,' density matrix elements between invariant',
     $            ' orbital sets have been zeroed ')
      end if
c
      if (debug) then
c        write(iout,*)' outgoing density matrix '
c        call print(den,naot,nao,iout)
      end if
c
      return
      end
