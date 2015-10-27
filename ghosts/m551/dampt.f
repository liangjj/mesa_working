*deck %W%  %G%
      subroutine dampt(ahess,smlld,shftd,nmx)
c
c***begin prologue     dampt
c***date written       871022   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c
c***keywords
c***author             lengsfield, byron (brl)
c***source             %W%   %G%
c
c***purpose            damp diagonal elements of the orbital hessian
c
c***description
c
c***references
c
c***routines called    (none)
c
c***end prologue       dampt
c
      implicit real*8 (a-h,o-z)
      dimension ahess(*)
c
      common /io/ inp,iout
c
      ix=1
      do 1 i=2,nmx
         if(ahess(ix+i).le.smlld) then
            damp=ahess(ix+i)+shftd
            write(iout,257) i,ahess(ix+i),damp
 257        format(/,'  damping orbital hessian row = ',i5,/,
     $               '  old diagonal ',f12.8,'   new diagonal ',f12.8)
            ahess(ix+i)=damp
         end if
 256     continue
         ix=ix+i
    1 continue
c
c
      return
      end
