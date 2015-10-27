*deck @(#)pm618.f	5.1  11/6/94
      subroutine pm618(zc,ac)
c***begin prologue     pm618.f
c***date written       yymmdd
c***revision date      11/6/94
c   february 10, 1994  rlm at lanl
c      isolating the SOBOL routines to generate the bsj surface grid.
c   january 7, 1994    rlm at lanl
c      mesa version 1.1 from version 6.0 compliments of gjt.
c***keywords           boundary surface, jello, solvent
c***author             tawa, greg and martin, richard(lanl)
c***source             @(#)pm618.f	5.1   11/6/94
c***purpose            generates the boundary surface for evaluating
c                      the electrostatic potential induced by a
c                      dielectric on a solute enclosed in a cavity.
c***description
c
c       This code generates the boundary surface for the BSJ calculation
c       of electrostatic potentials for the jello models
c
c       The code needs as input
c          (i) the centers and radii of a set of spheres,
c       the union of which represents the van der Waals surface of the
c       solute molecule,
c          (ii) the number of focus points on the solute van der Waals
c       surface where one desires to calculate a surface charge 
c       the surface charges at the location of the focus points 
c       represents the polarization of the solvent due to the charge
c       distribution of the solute
c           (iii) number of points to place on the solute van der Waals 
c       surface.
c
c***references
c
c***routines called
c
c***end prologue       pm618.f
      implicit none
c     --- input variables -----
c     --- input arrays (unmodified) ---
c     --- input arrays (scratch) ---
      real*8 zc(*)
      integer ac(*)
c     --- output arrays ---
c     --- output variables ---
c     --- scratch arrays ---
c     --- local variables ---
      character*4096 ops
      integer nat,intkey
      integer ncnt,npoints,nfocmx
      integer mpoints,mpoint2
      integer xcnt,ycnt,zcnt,rad
      integer c,ian,top
      integer wpadti,iadtwp,canget
      integer mark,nark
      integer xq,yq,zq,indx,u
      integer xyzch,ifocus
      integer icindx,ipindx,index1,index2,index3
      integer inp,iout
      integer kount,nfocus
      integer maxatm,i
      real*8 toang
      logical prnt,logkey
c
      common/io/inp,iout
c
      parameter(maxatm=2000,mpoints =40000)
      parameter(mpoint2 = 512)
      character*16 chdum(maxatm)
c
c
 1000 format(1x,'m618:generate cavity surface')
 1010 format(5x,'cavity radii(au)')
 1020 format(8x,a8,f13.4)
c
c
      write(iout,1000)
c
c     --- recover the options string ---
      call iosys('read character options from rwf',-1,0,0,ops)
c
c     --- get some basic information ---
      call iosys('read integer "number of atoms" from rwf',
     $           1,nat,0,' ')
c
c     --- check on options and set defaults ---
      npoints=intkey(ops,'solvent=npoints',mpoints,' ')
      ncnt=intkey(ops,'solvent=ncenters',nat,' ')
      nfocmx=intkey(ops,'solvent=nfocus',mpoint2,' ')
      prnt=logkey(ops,'solvent=print',.false.,' ')
c
c     --- allocate some core
      xcnt=1
      ycnt=xcnt+ncnt
      zcnt=ycnt+ncnt
      rad=zcnt+ncnt
      c=rad+ncnt
      ian=wpadti(c+3*nat)
      top=ian+nat
c
c     --- get atomic coordinates and atomic numbers ---
      call iosys('read real coordinates from rwf',-1,zc(c),0,' ')
      call iosys('read integer "atomic numbers" from rwf',
     $           nat,ac(ian),0,' ')
c
c     -- get some conversion constants ---
      call iosys('read real angstrom/bohr from rwf',1,toang,0,' ')
c
c     --- assign remaining input variables ---
c
c xcnt(i),ycnt(i),zcnt(i) = Coordinates of van der Waals sphere center i
c rad(i)                  = Radius of van der Waals sphere center i.
c       npoints           = Total number of points used to represent the
c                           molecular surface.
c       ncnt              = The number of van der Waals centers.
c       prnt              = Print flag.
c       nfocmx            = The number of non-buried  points to be used 
c                           as focus points.  If nfocmx = the number of
c                           non-buried points no averaging is performed.  
c                           If nfocmx < the total number of non-buried points 
c                           the code assigns all non-focus points to one focus
c                           point. This set of points makes up a plaquette and 
c                           the various matrix elements for each focus
c                           point will consist of an average over all of
c                           the points in the plaquette associated with the 
c                           focus point.
c
c     --- read the solvent radii.  if an element is different from zero,
c         it was read in the coordinate section and will not
c         be replaced by getrad. the others will get default values.
      call iosys('read real "solvent radii" from rwf',-1,zc(rad),0,' ')
      call getrad(zc(c),ac(ian),nat,zc(xcnt),zc(ycnt),zc(zcnt),zc(rad),
     $           ncnt,prnt,toang)
c
c     --- write the solvent radii to the rwf.
      call iosys('write real "solvent radii" to rwf',
     $            nat,zc(rad),0,' ')
c
c     --- print the radii ---
      call iosys('read character "z-names w/o dummies" from rwf',
     $           -1,0,0,chdum(1))
      write(iout,1010)
      do 10 i=1,nat
         write(iout,1020) chdum(i),zc(rad+i-1)
   10 continue
c
c     --- allocate more memory ---
      mark=ian
      nark=mark+ncnt+1
      xq=iadtwp(nark+ncnt+1)
      yq=xq+npoints
      zq=yq+npoints
      indx=wpadti(zq+npoints)
      u=iadtwp(indx+npoints)
c     --- arrays related to calculation of surface charges
      xyzch=u+2*nfocmx
      ifocus=wpadti(xyzch+3*nfocmx)
c     --- arrays related to clump algorithm
      icindx=ifocus+ncnt
      ipindx=icindx+npoints
      index1=ipindx+npoints
      index2=index1+nfocmx
      index3=index2+npoints
      top=index3+nfocmx
      call getscm(0,zc(1),canget,'m618',0)
      if(top.gt.canget) then
         write(iout,*) 'need',top,' integer words; have',canget
         call lnkerr('need too much core in m618')
      endif

      call surface(zc(xyzch),zc(xcnt),zc(ycnt),zc(zcnt),zc(rad),
     $         npoints,ncnt,
     $         prnt,nfocmx,ac(mark),ac(nark),
     $         zc(xq),zc(yq),zc(zq),ac(indx),zc(u),kount,
     $         nfocus,ac(ifocus),ac(icindx),ac(ipindx),ac(index1),
     $         ac(index2),ac(index3))
c
c     --- store surface coordinates ---
      call iosys('write integer "number of solvent surface points"'
     $           //' to rwf',1,nfocus,0,' ')
      call iosys('write real "solvent surface coordinates" to rwf',
     $           3*nfocus,zc(xyzch),0,' ')
c
c     --- exit gracefully ---
      call chainx(0)
c
c
      return
      end
