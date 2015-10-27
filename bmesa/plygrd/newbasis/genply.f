*deck region.f
c***begin prologue     region
c***date written       951229   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           
c***author             schneider, barry (nsf)
c***source             
c***purpose            lobatto polynomials.
c***                   
c***references         
c
c***routines called    
c***end prologue       region
      subroutine region(ppt,ppoly,ploc,ptq,ptp,xl,xr,gr2use,ngr,
     1                  ngrid,npt,npnts)
      implicit integer (a-z)
      real*8 pt, poly, ptreg, plyreg, xl, xr
      dimension gr2use(ngr), npt(ngrid), ploc(2,ngrid), npnts(ngr)
      common/io/inp, iout
      pointer (pptwt,pt(1))
      pointer (ppoly,poly(1))
      pointer (ptq,ptreg(1))
      pointer (ptp,plyreg(1))
      cnti=1
      cntj=1
      do 10 grd1=1,ngr
         grdi=g2use(grd1)
         npnts(grd1)=npt(grdi)
         q=cnti
         wt=q+npt(grdi)
         cnti=wt+npt(grdi)
         do 20 grd2=1,ngr
            grdj=g2use(grd2)
            p=cntj
            dp=p+npt(grd2)*grd(1)
            ddp=dp+npt(grd2)*grd(1)
            cntj=ddp+npt(grd2)*grd(1)
 20      continue
 10   continue   
         begin=ploc(1,grdi)
         
      cntpti=1
      do 10 grdi=1,ngrid
         qi=cntpti
         wti=qi+npt(grdi)
         cntpti=wti+npt(grdi)
         cntptj=1
         do 20 grdj=1,ngrid
            qj=cntptj
            wtj=qj+npt(grdj)
            cntpti=wtm+npt(grdj)
            p=cntply
            dp=p+npt(grdi)*npt(grdj)
            ddp=dp+npt(grdi)*npt(grdj)
            cntply=ddp+npt(grdi)*npt(grdj)
            call lgngr(poly(p),poly(dp),poly(ddp),pt(qi),pt(qj),
     1                 npt(grdi),npt(grdj),.false.)
 20      continue
 10   continue   
      return
      end       



