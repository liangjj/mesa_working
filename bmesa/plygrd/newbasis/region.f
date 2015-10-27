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
      subroutine region(px,pxply,py,pyply,xloc,pxloc,yloc,pyloc,
     1                  xl,xr,gr2use,ngr,ngrid,maxgrd,
     2                  npt,npnts,ngot)
      implicit integer (a-z)
      real*8 x, xply, y, yply, xl, xr, a, b, ainv, ainvsq
      dimension yloc(maxgrd), pyloc(maxgrd,maxgrd)
      dimension xloc(ngrid), pxloc(ngrid,ngrid)
      dimension gr2use(ngr), npt(ngrid), npnts(ngr)
      dimension ngot(2)
      common/io/inp, iout
      pointer (px,x(1))
      pointer (pxply,xply(1))
      pointer (py,y(1))
      pointer (pyply,yply(1))
      cnti=1
      cntji=1
      do 10 grdi=1,ngr
         grdxi=g2use(grdi)
         npnts(grdi)=npt(grdxi)
         yloc(grdi)=cnti
         yi=cnti
         ywti=yi+npt(grdxi)
         cnti=ywti+npt(grdxi)
         do 20 grdj=1,ngr
            grdxj=g2use(grdj)
            pyloc(j,i)=cntji
            pyji=cntji
            dpyji=pyji+npt(grdxi)*npt(grdxj)
            ddpyji=dpyji+npt(grdxi)*npt(grdxj)
            cntji=ddpyji+npt(grdxi)*npt(grdxj)
 20      continue
 10   continue   
      cnti=wpadti(cnti)
      cntji=wpadti(cntji)
      call memory(cnti,py,ngot(1),'y',0)
      call memory(cntji,pyply,ngot(2),'py',0)
      a=(xr-xl)*.5d0
      b=(xl+xr)*.5d0
      ainv=1.d0/a
      ainvsq=ainv*ainv
      do 30 grdi=1,ngr
         yi=yloc(grdi)
         ywti=yi+npnts(grdi)
         grdxi=g2use(grdi)
         xi=xloc(grdxi)
         xwti=xi+npt(grdxi)
         call shiftxy(x(xi),x(wti),y(yi),y(ywti),a,b,npnts(grdi))
 30   continue
      do 40 grdi=1,ngr
         grdxi=g2use(grdi)
         do 50 grdj=1,ngr
            grdxj=g2use(grdj)
            pxji=pxloc(grdxj,grdxi)   
            pyji=pyloc(grdj,grdi)   
            call copy(xply(pxji),yply(pyji),npnts(grdi)*npnts(grdj))
            dpxji=pxji+npnts(grdi)*npnts(grdj)
            dpyji=pyji+npnts(grdi)*npnts(grdj) 
            call vscale(yply(dpyji),xply(dpxji),ainv,
     1                  npnts(grdi)*npnts(grdj))
            ddpxji=dpxji+npnts(grdi)*npnts(grdj)
            ddpyji=dpyji+npnts(grdi)*npnts(grdj)
            call vscale(yply(ddpyji),xply(ddpxji),ainvsq,
     1                  npnts(grdi)*npnts(grdj)) 
 50      continue
 40   continue   
      return
      end       





