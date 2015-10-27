      subroutine printi(ipt,len,ii)
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
      dimension ipt(*)
c
      common /io/ inp,iout
c
c
      ix=0
      do 1 i=1,len
         write(iout,2)(ipt(ix+j),j=1,i)
         ix=ix+i
 1    continue
 2    format(5(2x,i8))
c
      return
      end
