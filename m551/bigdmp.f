*deck @(#)bigdmp.f	1.1  11/30/90
      subroutine bigdmp(ah,nmx,js,lb,smlld,shftd,diag)
c
c***begin prologue     bigdmp
c***date written       871022   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c
c***keywords
c***author             lengsfield, byron (brl)
c***source             @(#)bigdmp.f	1.1   11/30/90
c
c***purpose
c
c***description
c
c***references
c
c***routines called    (none)
c
c***end prologue       bigdmp
c
      implicit real*8 (a-h,o-z)
      dimension ah(nmx,lb),diag(lb)
c
      common /io/ inp,iout
c
c
      do 1 i=1,lb
c
         if(ah(js+i,i).gt.smlld)go to 256
         damp=ah(js+i,i)+shftd
         write(iout,257) i,ah(js+i,i),damp
 257     format(/,'  damping orbital hessian row = ',i5,/,
     $        '  old diagonal ',f12.8,'   new diagonal ',f12.8)
         ah(js+i,i)=damp
 256     continue
c
         diag(i)=ah(js+i,i)
 1    continue
c
      return
      end
