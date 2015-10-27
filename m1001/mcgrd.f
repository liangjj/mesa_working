*deck @(#)mcgrd.f	1.1  11/30/90
      subroutine mcgrd(grad,xlag,nbf,nco,nao,nvo,lmix,lmixt,
     $     mixhes,nob)
c
c***begin prologue
c***date written       871022   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c
c***keywords
c***author             lengsfield, byron (brl)
c***source             @(#)mcgrd.f	1.1   11/30/90
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
c     implicit real*8(a-h,o-p,r-z),             integer*2(q)
cc
cmp   extended dummy grad,xlag,mixhes
cc
      implicit real*8(a-h,o-z)
      dimension grad(2),xlag(nbf,2),mixhes(nob,2)
c
      common /io/ inp,iout
c
c
c
      lmix=0
      nocc=nco+nao
      nob=nocc+nvo
c
c
c     write(iout,8801)((xlag(j,i),j=1,nob),i=1,nocc)
c8801 format(' xlag ',/,5(2x,f16.8))
      if(nco.eq.0)go to 31
      do 30 i=1,nco
         if(nao.eq.0) go to 11
         js=nco+1
         do 10 j=js,nocc
            lmix=lmix+1
            grad(lmix)=xlag(j,i)-xlag(i,j)
 10      continue
 11      continue
         if(nvo.eq.0) go to 21
         js=nocc+1
         do 20 j=js,nob
            lmix=lmix+1
            grad(lmix)=xlag(j,i)
 20      continue
 21      continue
 30   continue
 31   continue
      if(nao.eq.0)go to 61
      is=nco+1
      do 60 i=is,nocc
         if(i.eq.nocc)go to 41
         js=i+1
         do 40 j=js,nocc
            if(mixhes(j,i).eq.0)go to 40
            lmix=lmix+1
            grad(lmix)=xlag(j,i)-xlag(i,j)
 40      continue
 41      continue
         if(nvo.eq.0) go to 51
         js=nocc+1
         do 50 j=js,nob
            lmix=lmix+1
            grad(lmix)=xlag(j,i)
 50      continue
 51      continue
 60   continue
 61   continue
c
      lmixt=lmixt+lmix
c
      return
      end
