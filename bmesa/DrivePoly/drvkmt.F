*deck drvkmt.f 
c***begin prologue     drvkmt
c***date written       000702   (yymmdd)
c***revision date               (yymmdd)
c***keywords           dvr
c***                   
c***author             schneider, b. i.(nsf)
c***source             drvkmt
c***purpose            driver for kinetic energy matrix elements
c
c***references       
c
c***routines called    iosys, util and mdutil
c***end prologue       drvkmt
      subroutine drvkmt(pt,wpt,cp,dcp,ddcp,kmat,cordsys,coord,
     1                  mass,edge,parity,n)
c
      implicit integer (a-z)
      character*(*) cordsys, coord, parity
      character*80 title
      real*8 pt, wpt, cp, dcp, ddcp, kmat, mass, edge, half
      logical prn
      dimension pt(n), wpt(n), cp(n,n), dcp(n,n), ddcp(n,n), kmat(n,n)
      dimension edge(2)
      common/io/inp, iout
      data prn/.true./
      data half / .5d0 /
c
      call kemel(kmat,cp,dcp,ddcp,pt,wpt,edge,n,cordsys,coord,parity)
      call vscale(kmat,kmat,-half/mass,n*n)
      if(prn) then
         title='kinetic energy matrix'
         call prntfm(title,kmat,n,n,n,n,iout)   
      endif
      return
      end


















