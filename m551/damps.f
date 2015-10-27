*deck @(#)damps.f	1.1  11/30/90
      subroutine damps(ahess,smlld,shftd,nmx)
c
c***begin prologue     damps
c***date written       871022   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c
c***keywords
c***author             lengsfield, byron (brl)
c***source             @(#)damps.f	1.1   11/30/90
c
c***purpose            damp diagonal elements of the orbital hessian
c
c***description
c
c***references
c
c***routines called    (none)
c
c***end prologue       damps
c
      implicit real*8 (a-h,o-z)
      logical debug
      dimension ahess(nmx,nmx)
c
      data debug/.false./
      common /io/ inp,iout
c
      do 1 i=2,nmx
         if(ahess(i,i).le.smlld) then
            damp=ahess(i,i)+shftd
            write(iout,257) i,ahess(i,i),damp
 257        format(/,'  damping orbital hessian row = ',i5,/,
     $               '  old diagonal ',f12.8,'   new diagonal ',f12.8)
            ahess(i,i)=damp
         end if
    1 continue
c
c
      if (debug) then
         write (iout,*) '  damps: ahess'
         call matout(ahess,nmx,nmx,nmx,nmx,iout)
      end if
c
      return
      end
