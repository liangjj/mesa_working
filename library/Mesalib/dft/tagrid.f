*deck @(#)tagrid.f	5.2 2/5/95
      subroutine tagrid(xyzgrid,grdwts,atnum,atcen,maxr,mxgrd,
     $                  nr,ngrid,ptrad,grdtyp,tmprpts,tmprwts,
     $                  mxang,tmppts,tmpwts)
c***begin prologue     tagrid.f
c***date written       950120  
c***revision date      2/5/95      
c***keywords           numerical integration, standard grid
c***author             martin, richard (lanl)
c***source             @(#)tagrid.f	5.2   2/5/95
c***purpose            returns the points and weights associated with
c                      the standard grids pruned according to
c                      treutler and alhrichs.
c***description
c
c***references
c                     O. Treutler and R. Ahlrichs, J.Chem. Phys., 346 (1995).
c
c***end prologue      tagrid.f
      implicit none
c
c     --- input variables -----
      integer atnum,maxr,mxgrd,mxang
      character*(*) grdtyp
c     --- input arrays (unmodified) ---
      real*8 atcen(3)
c     --- input arrays (scratch) ---
      real*8 tmprwts(maxr),tmprpts(maxr)
      real*8 tmpwts(mxang),tmppts(mxang,3)
c     --- output arrays ---
      integer ptrad(0:maxr)
      real*8 xyzgrid(mxgrd,3),grdwts(mxgrd)
c     --- output variables ---
      integer nr,ngrid
c     --- scratch arrays ---
c     --- local variables ---
      integer angsiz
      integer maxw,nang
      integer i,j,lebord
      integer inp,iout
      real*8 xi(104)
      real*8 cut1,cut2,one,two,three
      real*8 scale,toang
      logical called
c
      parameter (one=1.0d+00,two=2.0d+00,three=3.0d+00)
      data xi/
c                       H    He
     $              0.8d0,0.9d0,
c                      Li    Be     B     C     N     O     F    Ne
     $              1.8d0,1.4d0,1.3d0,1.1d0,0.9d0,0.9d0,0.9d0,0.9d0,
c                      Na    Mg    Al    Si     P     S    Cl    Ar
     $              1.4d0,1.3d0,1.3d0,1.2d0,1.1d0,1.0d0,1.0d0,1.0d0,
c                       K    Ca
     $              1.5d0,1.4d0,
c                      Sc    Ti     V    Cr    Mn
     $              1.3d0,1.2d0,1.2d0,1.2d0,1.2d0,
c                                  Fe    Co    Ni    Cu    Zn
     $                          1.2d0,1.2d0,1.1d0,1.1d0,1.1d0,
c                                  Ga    Ge    As    Se    Br    Kr
     $                          1.1d0,1.0d0,0.9d0,0.9d0,0.9d0,0.9d0,
     $           68*0.0d0/
      data called/.false./
c
      save xi,toang,called
      common /io/inp,iout
c
c
      if(.not.called) then
         call iosys('read real angstrom/bohr from rwf',1,toang,0,' ')
         called=.true.
      endif
c
c     --- check to see if the grids are available.
      if(atnum.gt.36) then
         call lnkerr(
     $        'Ack.  You called tagrid for something above krypton')
      endif
c
c     --- set up basic quadrature parameters depending on grid and element.
      if(grdtyp.eq.'tag1') then
         if (atnum.le. 2) then
c           H - He
            nr=20
            maxw=11
         else if (atnum .le. 10) then
c           Li-Ne
            nr=25
            maxw=17
         else if (atnum .le. 18) then
c           Na-Ar
            nr=30
            maxw=17
         else if(atnum.le.36) then
c           K-Kr
            nr=35
            maxw=17
         endif
      else if(grdtyp.eq.'tag2') then
         if (atnum.le. 2) then
c           H - He
            nr=25
            maxw=17
         else if (atnum .le. 10) then
c           Li-Ne
            nr=30
            maxw=23
         else if (atnum .le. 18) then
c           Na-Ar
            nr=35
            maxw=23
         else if(atnum.le.36) then
c           K-Kr
            nr=40
            maxw=23
         endif
      else if(grdtyp.eq.'tag3') then
         if (atnum.le. 2) then
c           H - He
            nr=30
            maxw=23
         else if (atnum .le. 10) then
c           Li-Ne
            nr=35
            maxw=29
         else if (atnum .le. 18) then
c           Na-Ar
            nr=40
            maxw=29
         else if(atnum.le.36) then
c           K-Kr
            nr=45
            maxw=29
         endif
      else if(grdtyp.eq.'tag4') then
         if (atnum.le. 2) then
c           H - He
            nr=35
            maxw=29
         else if (atnum .le. 10) then
c           Li-Ne
            nr=40
            maxw=35
         else if (atnum .le. 18) then
c           Na-Ar
            nr=45
            maxw=35
         else if(atnum.le.36) then
c           K-Kr
            nr=50
            maxw=35
         endif
      else if(grdtyp.eq.'tag5') then
         if (atnum.le. 2) then
c           H - He
            nr=45
            maxw=35
         else if (atnum .le. 10) then
c           Li-Ne
            nr=50
            maxw=47
         else if (atnum .le. 18) then
c           Na-Ar
            nr=55
            maxw=47
         else if(atnum.le.36) then
c           K-Kr
            nr=60
            maxw=47
         endif
      else
         write(iout,*) grdtyp
         call lnkerr('unrecognized grid to tagrid')
      endif
c
c     --- generate radial quadrature
      ngrid=0
      scale=xi(atnum)/toang
      call alrixq(tmprpts,tmprwts,nr,scale)
c
c 
c     --- prune the grid.  use l=5 angular quadrature for points near
c         the nucleus, l=11 for intermediate, and maxw for the rest.
c         the quadrature generation puts large r at the beginning,
c         so reverse the loop here to put points near r=0 first.
      cut1=float(nr)/three
      cut2=float(nr)/two
      ptrad(0)=1
      do 10 i=nr,1,-1
         if(float(nr-i+1).le.cut1) then
c           use the icosahedral quadrature for l=5
            lebord=5
            nang=12
         else if(float(nr-i+1).le.cut2) then
            lebord=11
            nang=angsiz(lebord)
         else
            lebord=maxw
            nang=angsiz(lebord)
         endif
         call sphere(mxang,tmppts,tmpwts,lebord,nang)
         do 20 j=1,nang
            xyzgrid(j+ngrid,1)=tmppts(j,1)*tmprpts(i)+atcen(1)
            xyzgrid(j+ngrid,2)=tmppts(j,2)*tmprpts(i)+atcen(2)
            xyzgrid(j+ngrid,3)=tmppts(j,3)*tmprpts(i)+atcen(3)
            grdwts(j+ngrid)=tmpwts(j)*tmprwts(i)
 20      continue 
         ptrad(i)=ptrad(i-1)+nang
         ngrid=ngrid+nang
 10   continue 
c
c
      if (ngrid .gt. mxgrd) then
         write(iout,*) 'tagrid: ngrid exceeds mxgrid; ngrid,mxgrid',
     $                  ngrid,mxgrd
         write(iout,*) 'resubmit with scf=(mxgrid=',ngrid,')'
         call lnkerr('tagrid: ngrid exceeds mxgrd.')
      endif
c
c
      return
      end

