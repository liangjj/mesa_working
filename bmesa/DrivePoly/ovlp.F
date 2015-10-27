*deck ovlp.f 
c***begin prologue     ovlp
c***date written       000702   (yymmdd)
c***revision date               (yymmdd)
c***keywords           dvr
c***                   
c***author             schneider, b. i.(nsf)
c***source             ovlp
c***purpose            overlap matrix
c
c***references       
c
c***routines called    iosys, util and mdutil
c***end prologue       ovlp
      subroutine ovlp(pntxwt,n,typwt,typarg)
c
      implicit integer (a-z)
      character*(*) typwt, typarg
      real*8 xwt
      common/io/inp, iout
      pointer (pntxwt,xwt(1))
c
      pt=1
      wpt=pt+n
      cp=wpt+n
      dcp=cp+n*n
      ddcp=dcp+n*n
      call renorm(xwt(cp),xwt(dcp),xwt(ddcp),xwt(wpt),n)
      return
      end


















