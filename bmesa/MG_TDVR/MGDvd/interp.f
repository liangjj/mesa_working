*deck interp.f 
c***begin prologue     interp
c***date written       960718   (yymmdd)
c***revision date               (yymmdd)
c***keywords           time, dvr, orthogonal polynomial
c***                   bose-einstein
c***author             schneider, b. i.(nsf)
c***source             tprop
c***purpose            interpolation of a vector from gridi to gridj 
c***                   using the DVR functions. 
c***                   
c***description        
c***references       
c
c***routines called    iosys, util and mdutil
c***end prologue       mgdrv
      subroutine interp(pm,psln,pbasis,psoln,nmax,ntot,m,grdi,grdj,
     1                  mgr,nvc,dim)
c
      implicit integer (a-z)
      real*8 basis, soln
      dimension pm(mgr,mgr,dim), psln(mgr), nmax(4,mgr), ntot(mgr)
      dimension m(mgr)
      common/io/inp, iout
      pointer(pbasis,basis(1))
      pointer(psoln,soln(1))
      soli=psln(grdi)         
      nxi=nmax(2,grdi)
      nyi=nmax(3,grdi)
      nzi=nmax(4,grdi)
      nti=nmax(1,grdi)
      ni=ntot(grdi)
      mi=m(grdi)
      pxii=pm(grdi,grdi,2)
      pyii=pm(grdi,grdi,3)
      pzii=pm(grdi,grdi,4)
      ptii=pm(grdi,grdi,1)
      solj=psln(grdj)
      nxj=nmax(2,grdj)
      nyj=nmax(3,grdj)
      nzj=nmax(4,grdj)
      ntj=nmax(1,grdj)
      nj=ntot(grdj)
      mj=m(grdj)
      pxij=pm(grdi,grdj,2)
      pyij=pm(grdi,grdj,3)
      pzij=pm(grdi,grdj,4)
      ptij=pm(grdi,grdj,1)
      call gi2gj(soln(soli),soln(solj),
     1           basis(ptij),basis(pxij),basis(pyij),basis(pzij),
     2           basis(ptii),basis(pxii),basis(pyii),basis(pzii),
     3           nti,nxi,nyi,nzi,ntj,nxj,nyj,nzj,ni,nj,nvc,dim)
      call vivo(soln(soli),soln(soli),
     1          basis(ptii),basis(pxii),basis(pyii),basis(pzii),
     2          nti,nxi,nyi,nzi,ni,nvc,dim)
      return
      end







