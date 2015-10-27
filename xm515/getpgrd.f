*deck %W% %G%
      subroutine getpgrd( )
c***begin prologue     %M%
c***date written       951121      (yymmdd)  
c***revision date      %G%
c
c***keywords           
c***author             russo, thomas(lanl)
c***source             %W% %G%
c***purpose            generate a grid for the poisson solver
c                      
c***description        
c
c     general grid is simple: just define one radshl block with all of the
c           radial points in it, all of them having the same lebord.
c     "sg1" grid is like an sg1 grid for integration, but uses cotes points
c           and weights instead of euler maclaurin.  So we need a new routine,
c           since sg1 calls eulerq to get radial grid.  Icky-tuu, but
c           c'est la vie.  The upside is that psg1 can eventually be used
c           to replace the sg1 code, but that'll mean that some code will
c           need to be rewritten.
c     Other grids will be forthcoming, if someone else codes them.
c
c***references
c
c***routines called
c
c***end prologue       %M%
c
      implicit none
c --- input variables
c --- input arrays (unmodified)
c --- input arrays (changed)
c --- output variables
c --- output arrays
c --- scratch arrays
c --- local storage
      character*32 grdnm
 900  format('atom ',i3.3,' rgrd')
 910  format('atom ',i3.3,' rwt')
 1000 format('atom ',i3.3,' nrblk')
 1010 format('atom ',i3.3,' rblkord')
 1020 format('atom ',i3.3,' rblkbgn')
 1025 format('atom ',i3.3,' rblksiz')
 1030 format('atom ',i3.3,' agrd ',i3.3)
 1040 format('atom ',i3.3,' angwt ',i3.3)
 1050 format('atom ',i3.3,' nlm')
 1050 format('atom ',i3.3,' nlmx')
 1060 format('atom ',i3.3,' ylm ',i3.3)
 1070 format('atom ',i3.3,' pgrid')
 1080 format('atom ',i3.3,' pwts')
 1089 format('atom ',i3.3,' vwts')
      ylmermx=0.0
      do 100 iatom=1,natoms
         if (pgrdtyp(iatom).eq.'general') then
            nrings=(vradial/(vncrule-1))+1
            rhomx=rhomax('bragg',ian(iatom))
c           lay down radial grid and weights & save
            call cotes(rhomx,vradial,vncrule,nrings,rpts,rwt,jacob,
     $           npring,nrings,nr,nwtot)
            write(unit=grdnm,fmt=900)iatom
            call iosys('write real "'//grdnm//'" to '//grdfil,
     $           nr,rpts,0,' ')
            write(unit=grdnm,fmt=910)iatom
            call iosys('write real "'//grdnm//'" to '//grdfil,
     $           nr,rwt,0,' ')
c           lay down angular quadrature
            call sphere(vnomega,angpts,angwt,lebord,vnomega)
c           save angular grids (trivial entries in data structure)
            nrblk=1
            rblkord(1)=lebord
            rblkbgn(1)=1
            rblksiz(1)=nr
            write(unit=grdnm,fmt=1000)iatom
            call iosys('write integer "'//grdnm//'" to '//grdfil,
     $           1,nrblk,0,' ')
            write(unit=grdnm,fmt=1010)iatom
            call iosys('write integer "'//grdnm//'" to '//grdfil,
     $           1,rblkord,0,' ')
            write(unit=grdnm,fmt=1020)iatom
            call iosys('write integer "'//grdnm//'" to '//grdfil,
     $           1,rblkbgn,0,' ')
            write(unit=grdnm,fmt=1025)iatom
            call iosys('write integer "'//grdnm//'" to '//grdfil,
     $           1,rblksiz,0,' ')
            write(unit=grdnm,fmt=1030)iatom,1
            call iosys('write real "'//grdnm//'" to '//grdfil,
     $           vnomega*3,angpts,0,' ')
            write(unit=grdnm,fmt=1040)iatom,1
            call iosys('write real "'//grdnm//'" to '//grdfil,
     $           vnomega,angwt,0,' ')
c           form ylms on angular grid & save, test orthonormality to get
c           ylmerr
            nlm(1)=(vlmax+1)*(vlmax+1)
            write(unit=grdnm,fmt=1050)iatom
            call iosys('write integer "'//grdnm//'" to '//grdfil,
     $           1,nlm,0,' ')
c for this case nlm and nlmx are same
            write(unit=grdnm,fmt=1055)iatom
            call iosys('write integer "'//grdnm//'" to '//grdfil,
     $           1,nlm,0,' ')
            call vylm(vlmax,vnomega,0,angpts,test,vnomega,nlm(1),
     $           ylm,ptlm,ctheta,phi,cphi,sphi,plm,scr)
            call tstylm(vlmax,vnomega,nlm(1),ylm,angwt,ptlm,scr,
     $           ylmerr)
            ylmermx=max(ylmermx,ylmerr)
            write(unit=grdnm,fmt=1060)iatom,1
            call iosys('write real "'//grdnm//'" to '//grdfil,
     $           nlm(1)*vnomega,ylm,0,' ')
c           form catenated grid&save with weights
            call catenat(nr,ngrid,vnomega,vnomega,rpts,rwt,
     $           angpts,angwt,grid,gridwts,vptrad,c(1,iatom))
            write(unit=grdnm,fmt=1070)iatom
            call iosys('write real "'//grdnm//'" to '//grdfil,
     $           ngrid*3,grid,0,' ')
            write(unit=grdnm,fmt=1080)iatom
            call iosys('write real "'//grdnm//'" to '//grdfil,
     $           ngrid,gridwts,0,' ')
         else if (pgrdtyp(iatom).eq.'sg1') then
            radial=51
            ncrule=5
            nrings=(radial/(ncrule-1))+1
            rhomx=rhomax('bragg',ian(iatom))
c        ---lay down radial quadrature and save
            call cotes(rhomx,radial,ncrule,nrings,rpts,rwt,jacob,
     $           npring,nrings,nr,nwtot)
            write(unit=grdnm,fmt=900)iatom
            call iosys('write real "'//grdnm//'" to '//grdfil,
     $           nr,rpts,0,' ')
            write(unit=grdnm,fmt=910)iatom
            call iosys('write real "'//grdnm//'" to '//grdfil,
     $           nr,rwt,0,' ')
c        ---determine blocks for pruning grids:
C              psg1 takes radial grid and this atnum, returns:
c              nrblk: number of radial shell blocks on this atom (always 5
c                for the sg1 grids)
c              angmx: number of angular points on largest grid in this atom
c                   (always 23 for sg1)
c              nlmx: maximum #ylms on any grid for this atom
c                   (always (11+1)*(11+1) = 144 corresponding to leb=23)
c              rblkord(array of nrblk): lebedev order of each block
c              rblkbgn(array of nrblk): starting radial grid point this block
c              rblksiz(array of nrblk): number grid points this block
c              nlm(array of nrblk): number of ylms needed 
            call psg1(nr,rpts, ian(iatom), rhomx, nrblk,angmx, nlmx,
     $           rblkord, rblkbgn, rblksiz,nlm)
            write(unit=grdnm,fmt=1000)iatom
            call iosys('write integer "'//grdnm//'" to '//grdfil,
     $           1,nrblk,0,' ')
            write(unit=grdnm,fmt=1010)iatom
            call iosys('write integer "'//grdnm//'" to '//grdfil,
     $           nrblk,rblkord,0,' ')
            write(unit=grdnm,fmt=1020)iatom
            call iosys('write integer "'//grdnm//'" to '//grdfil,
     $           nrblk,rblkbgn,0,' ')
            write(unit=grdnm,fmt=1025)iatom
            call iosys('write integer "'//grdnm//'" to '//grdfil,
     $           nrblk,rblksiz,0,' ')
c           lay down each angular quadrature & save (data structure for
c                  retrieval)
            ngrid=0
            do 200 j=1,nrblk
               lebord=rblkord(j)
               nang=angsiz(lebord)
               call sphere(mxang,angpts,angwt,lebord,nang)
               write(unit=grdnm,fmt=1030)iatom,j
               call iosys('write real "'//grdnm//'" to '//grdfil,
     $              nang*3,angpts,0,' ')
               write(unit=grdnm,fmt=1040)iatom,j
               call iosys('write real "'//grdnm//'" to '//grdfil,
     $              nang,angwt,0,' ')
c           form ylms this each angular quadrature & save
               lmax=(lebord-1)/2
               nlm(j)=(lmax+1)**2
               call vylm(lmax,nang,ioff,angpts,test,nang,
     $              nlm(j),ylm,ptlm,ctheta,phi,cphi,sphi,plm,scr)
               call tstylm(lmax,nang,nlm(j),ylm,angwt,ptlm,scr,
     $              ylmerr)
               ylmermx=max(ylmerr,ylmermx)
               write(unit=grdnm,fmt=1060)iatom,j
               call iosys('write real "'//grdnm//'" to '//grdfil,
     $              nlm(j)*nang,ylm,0,' ')
c           form catenated grid&save with weights
               do 250 k=rblkbgn(j),rblkbgn(j)+rblksiz(j)-1
                  do 275 l=1,nang
                     grid(ngrid+l,1)=angpts(l,1)*rpts(k)+atcen(1,iatom)
                     grid(ngrid+l,2)=angpts(l,2)*rpts(k)+atcen(2,iatom)
                     grid(ngrid+l,3)=angpts(l,3)*rpts(k)+atcen(3,iatom)
                     gridwts(ngrid+l)=angwt(l)*rwt(k)
 275              continue 
                  ngrid=ngrid+nang
 250           continue 
               write(unit=grdnm,fmt=1070)iatom
               call iosys('write real "'//grdnm//'" to '//grdfil,
     $              ngrid*3,grid,0,' ')
               write(unit=grdnm,fmt=1080)iatom
               call iosys('write real "'//grdnm//'" to '//grdfil,
     $              ngrid,gridwts,0,' ')
 200        continue 
         else
            call plnkerr('Unknown grid type in getpgrd',5566)
         endif
c        form voronoi weights for catenated grid & save
         call voronoi(natoms,c,grid,wts,ngrid,ngrid,iatom,
     $        vwts,rnuc,amu,pwtx,rr,adjust,radii)
         write(unit=grdnm,fmt=1090) iatom
         call iosys('write real "'//grdnm//'" to '//grdfil,
     $        ngrid,vwts,0,' ')
 100  continue 
      return
      end
