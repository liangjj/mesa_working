      subroutine mtest(mix,mixo,isymm,lok,len,locsym,nocc,nsym)
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
cc
cmp   extended dummy mix,mixo,lok,len
cc
      dimension mix(2),mixo(2),isymm(2),lok(2),len(2),locsym(2),
     $     nocc(2)
c
      common /io/ inp,iout
c
c
      do 20 ii=1,nsym
         locii=locsym(ii)
         nocii=nocc(ii)
         isii=isymm(ii)
         write(iout,30) locii,nocii,isii
         do 10 jj=1,nocii
            locj=lok(locii+jj)
            lenj=len(locii+jj)
            write(iout,31)(mix(locj+isii+mm),mm=1,lenj)
            write(iout,32)(mixo(locj+isii+mm),mm=1,lenj)
 10      continue
 20   continue
 30   format(' locii  nocii  isii ',3i6)
 31   format('  mix  ',16i3)
 32   format('  mixo ',16i3)
      return
      end
