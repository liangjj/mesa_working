*deck @(#)mcgrad.f	5.1  11/6/94
      subroutine mcgrad(gradr,grads,ncor,ncos,naor,naos,nobr,nobs,
     $     cr,cs,nbfr,nbfs,hess,mr,ms,gcol,tempr,temps)
c
c***begin prologue
c***date written       871022   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c
c***keywords
c***author             lengsfield, byron (brl)
c***source             @(#)mcgrad.f	5.1   11/6/94
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
      implicit real*8(a-h,o-z)
cc
cmp   extended dummy gradr,grads,cr,cs,hess,gcol,tempr,temps
cc
      dimension gradr(nbfr,2),grads(nbfs,2),cr(nbfr,2),cs(nbfs,2),
     $     hess(2),gcol(2),tempr(nbfr,2),temps(nbfs,2)
c
      common / number / zero,pt5,one,two,four,eight
c
      mrs=nbfr*nbfs
      irs=ncor+1
      ire=ncor+naor
      iss=ncos+1
      ise=ncos+naos
c
      if(mr.ne.ms)go to 2000
c
      if(naor.eq.0)go to 4000
c
c

      lh=1
      do 1075 ir=irs,ire
c
         do 1070 jr=irs,ir
c
            call apbc(gradr(1,jr),hess(lh),cr(1,ir),nbfr,nbfr,1)
            call apbtc(gradr(1,ir),hess(lh),cr(1,jr),nbfr,nbfr,1)
            lh=lh+mrs
c
 1070    continue
 1075 continue
c
c
      go to 4000
c
 2000 continue
c
      if(naor.eq.0)go to 2096
c
 2060 continue
c

      lh=1
      if(naos.eq.0)go to 2076
      do 2075 is=iss,ise
c
         do 2070 jr=irs,ire
c
            call apbc(gradr(1,jr),hess(lh),cs(1,is),nbfr,nbfs,1)
            lh=lh+mrs
c
 2070    continue
 2075 continue
c
 2076 continue
c
cc
 2096 continue
cc
      if(naos.eq.0)go to 4000
cc
c

      lh=1
      if(naor.eq.0)go to 3076
      do 3075 is=iss,ise
c
         do 3070 jr=irs,ire
c
c
            call apbtc(grads(1,is),hess(lh),cr(1,jr),nbfr,nbfs,1)
            lh=lh+mrs
c
 3070    continue
 3075 continue
c
 3076 continue
c
c
c
 4000 continue
c
      return
      end
