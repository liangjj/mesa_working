*deck @(#)pm621.f	5.1  11/6/94
      subroutine pm621(zc,ac)
c***begin prologue     pm621.f
c***date written       yymmdd
c***revision date      11/6/94
c   january 7, 1994    rlm at lanl
c      mesa version 1.1 from version 6.0 compliments of gjt.
c***keywords           boundary surface, jello, solvent
c***author             tawa, greg and martin, richard(lanl)
c***source             @(#)pm621.f	5.1   11/6/94
c***purpose            evaluates electrostatic potential induced by a
c                      dielectric on a solute enclosed in a cavity.
c***description
c
c       This is the main program for the BSJ calculation of
c	electrostatic potentials for the jello models
c
c       The code needs as input (i) a set of source charges, one source 
c       charge for each solute atom, (ii) the centers and radii of a set of
c       spheres the union of which represents the van der Waals surface 
c       of the solute molecule, (iii) the locations of a set of observation 
c       points where the electrostatic potential due to the polarization 
c       of the solvent is to be calculated, (iv) the dielectric constant of 
c       the solvent, (v) number of focus points on the solute van der Waals
c       surface one desires to calculate a surface charge (the surface charges
c       at the location of the focus points represents the polarizatin of the 
c       solvent due to the charge distribution of the solute), 
c       and (vi) the total number of points to place on the solute vdW surface.
c
c
c       We save the x,y,z coordinates and the magnitudes of the point 
c       charges on the van der Waals surface of the solute on the rwf.
c       The x,y,z coordinates are those of the focus points.
c
c
c       The observation points are assumed to be interior to
c       the molecular volume where the dielectric susceptibility
c       is zero.
c
c       The dielectric constant is assumed to be a constant, bulk
c       value outside the molecular surface and one (1.0) inside.
c***references
c
c***routines called
c
c***end prologue       pm621.f
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
      character*3 answer
      integer nat,intkey
      integer nsrcs,ncnt,npoints,nobs,nfocus
      integer mpoints,mpoint2
      integer x,y,z,q,xcnt,ycnt,zcnt,rad,xobs,yobs,zobs
      integer c,ian,top
      integer mem,wpadti,iadtwp,canget
      integer mark,nark,phi
      integer ex,ey,ez,efgxx,efgxy,efgxz,efgyy,efgyz,efgzz
      integer xq,yq,zq,indx,u
      integer a2,b2,xyzch,charge,ifocus
      integer icindx,ipindx,index1,index2,index3
      integer ncharge
      integer inp,iout
      integer nbf
      integer qd,qq
      real*8 epsilon,fpkey
      real*8 aben,deltag2,pomf2,self
      real*8 toang,toev,tokcal,jph,jpcal,avog,etoesu,slight,fines,emass
      real*8 energ,senerg,thresh
      logical cnvg,prnt,logkey
      logical first
c
      common/io/inp,iout
c
      parameter(mpoints =40000)
      parameter(mpoint2 = 512)
c
c     --- default threshold for convergence in solvent field ---
      data thresh/1.0d-06/
c
 1000 format(1x,'m621:generate solvent potential')
 1010 format(5x,'solvent calculation converged')
 1020 format(5x,'gas phase energy          ',f15.9,
     $      /5x,'solvated energy           ',f15.9,
     $      /5x,'solvation energy          ',f15.9,
     $      /5x,'   classical:             ',
     $      /5x,'   solvation free energy  ',f15.9,
     $      /5x,'   potential of mean force',f15.9,
     $      /5x,'   solvent self energy    ',f15.9)
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
      call iosys('read integer "number of basis functions" from rwf',
     $           1,nbf,0,' ')
c
c     --- check on options and set defaults ---
      epsilon=fpkey(ops,'solvent=epsilon',epsilon,' ')
      npoints=intkey(ops,'solvent=npoints',mpoints,' ')
      nsrcs=intkey(ops,'solvent=nsources',nat,' ')
      ncnt=intkey(ops,'solvent=ncenters',nat,' ')
      nobs=intkey(ops,'solvent=nobs',nat,' ')
      nfocus=intkey(ops,'solvent=nfocus',mpoint2,' ')
      thresh=fpkey(ops,'solvent=convergence',thresh,' ')
      prnt=logkey(ops,'solvent=print',.false.,' ')
c
c     --- see if this is the first time through the iterative
c         solution in the solvent field.
      call iosys('read real "energy" from rwf',1,energ,0,' ')
      call iosys('does "solvated energy" exist on rwf',
     $            0,0,0,answer)
      if(answer.eq.'no') then
         first=.true.
c        --- then the current energy is the gas-phase energy
         aben=energ
         call iosys('write real "solute gas phase energy" to rwf',
     $              1,aben,0,' ')
      else
         first=.false.
c        --- test convergence ---
         cnvg=.false.
         call iosys('read real "solvated energy" from rwf',
     $              1,senerg,0,' ')
         call iosys('read real "solute gas phase energy"'
     $              //' from rwf',1,aben,0,' ')
         if(abs(senerg-energ).lt.thresh) then
            cnvg=.true.
         endif
      endif
      call iosys('write real "solvated energy" to rwf',
     $           1,energ,0,' ')
c
c     --- if converged, clean up some things
c         if not converged, then write the current scf solution
c            to the rwf as the guess vector ---
      if(cnvg) then
         write(iout,1010)
         call iosys('read real "solvent free energy" from rwf',
     $              1,deltag2,0,' ')
         call iosys('read real "solvent potential of mean force"'
     $              //' from rwf',1,pomf2,0,' ')
         call iosys('read real "solvent self energy" from rwf',
     $              1,self,0,' ')
         write(iout,1020) aben,energ,energ-aben,
     $                    deltag2,pomf2,self
c        NOTE THAT THIS IS AN EFFECTIVE EXIT
         goto 200
      else
         call iosys('read real "scf vector" from rwf',
     $              nbf*nbf,zc(1),0,' ')
         call iosys('write real "guess vector" to rwf',
     $              nbf*nbf,zc(1),0,' ')
      endif
c
c     --- allocate some core ---
      x=1
      y=x+nsrcs
      z=y+nsrcs
      q=z+nsrcs
      qd=q+nsrcs
      qq=qd+3*nsrcs
      xcnt=qq+6*nsrcs
      ycnt=xcnt+ncnt
      zcnt=ycnt+ncnt
      rad=zcnt+ncnt
      xobs=rad+ncnt
      yobs=xobs+nobs
      zobs=yobs+nobs
      c=zobs+nobs
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
      call iosys('read real j/hartree from rwf',1,jph,0,' ')
      call iosys('read real j/cal from rwf',1,jpcal,0,' ')
      call iosys('read real avogadro from rwf',1,avog,0,' ')
      call iosys('read real esu/e- from rwf',1,etoesu,0,' ')
      call iosys('read real lightspeed from rwf',1,slight,0,' ')
      call iosys('read real "electron mass" from rwf',1,emass,0,' ')
      call iosys('read real "fine-structure" from rwf',
     $            1,fines,0,' ')
      toev=(1.0d-05)*fines*fines*emass*slight*slight*slight/etoesu
      tokcal=jph*avog/(1.0d03*jpcal)
c
c     --- assign remaining input variables ---
c
c x(i),y(i),z(i)          = x,y,z coordinates of source charge i
c q(i)                    = Value of the source charge i.
c xcnt(i),ycnt(i),zcnt(i) = Coordinates of van der Waals sphere center i
c rad(i)                  = Radius of van der Waals sphere center i.
c xobs(i),yobs(i),zobs(i) = Coordinates of the observation point i. Note
c                           to calculate electrostatic contribution to the 
c                           free energy of solvation, the coordinates of the
c                           observation points must be identical to the
c                           coordinates of the source points.
c       epsilon           = Dielectric constant of the solvent.
c       npoints           = Total number of points used to represent the
c                           molecular surface.
c       nsrcs             = The number of sources.
c       ncnt              = The number of van der Waals centers.
c       nobs              = The number of observation points.
c       prnt              = Print flag.
c       aben              = Gas phase energy of the solute.
c       nfocus            = The number of non-buried  points to be used 
c                           as focus points.  If nfocus = the number of 
c                           non-buried points,no averaging is performed.
c                           If nfocus is less than the total number of 
c                           non-buried points the code assignes all non-focus
c                           points to one focus point.
c                           This set of points makes up a plaquette and 
c                           the various matrix elements for each focus
c                           point will consist of an average over all of
c                           the points in the plaquette associated with the 
c                           focus point.
      call input(zc(c),ac(ian),nat,zc(x),zc(y),zc(z),zc(q),
     $           zc(qd),zc(qq),zc(xcnt),zc(ycnt),zc(zcnt),zc(rad),
     $           zc(xobs),zc(yobs),zc(zobs),epsilon,npoints,nsrcs,
     $           ncnt,nobs,prnt,nfocus,toang,ops)
c
c     --- echo check input ---
      if(first) then
         call setup(zc(x),zc(y),zc(z),zc(q),zc(qd),zc(qq),
     $              zc(xcnt),zc(ycnt),zc(zcnt),zc(rad),
     $              zc(xobs),zc(yobs),zc(zobs),epsilon,npoints,nsrcs,
     $              ncnt,nobs,prnt,nfocus,thresh)
      endif
c
c     --- have we enough memory ---
c     mem = 4*msrcs
c    &    + 13*mobs
c    &    + 6*mcnt
c    &    + 3*mcntp1
c    &    + 7*mpoints +4*mpoint2+mpoint2*mpoint2
      call getscm(0,zc(1),canget,'m621',0)
c     mem = 4*nsrcs
c    &    + 13*nobs
c    &    + 6*ncnt
c    &    + 3*(ncnt+1)
c    &    + 7*npoints +4*nfocus+nfocus*nfocus
c     if(iadtwp(mem).gt.canget) then
c        write(iout,*) 'need',iadtwp(mem),' have',canget
c        call lnkerr('not enough memory to continue m621')
c     endif
c
c     --- allocate more memory ---
      mark=ian
      nark=mark+ncnt+1
      phi=iadtwp(nark+ncnt+1)
      ex=phi+nobs
      ey=ex+nobs
      ez=ey+nobs
      efgxx=ez+nobs
      efgxy=efgxx+nobs
      efgxz=efgxy+nobs
      efgyy=efgxz+nobs
      efgyz=efgyy+nobs
      efgzz=efgyz+nobs
c
      xq=efgzz+nobs
      yq=xq+npoints
      zq=yq+npoints
      indx=wpadti(zq+npoints)
      u=iadtwp(indx+npoints)
c     --- arrays related to calculation of surface charges
      a2=u+2*ncnt
      b2=a2+nfocus*nfocus
      charge=b2+nfocus
      xyzch=charge+nfocus
      ifocus=wpadti(xyzch+3*nfocus)
c     --- arrays related to clump algorithm
      icindx=ifocus+ncnt
      ipindx=icindx+npoints
      index1=ipindx+npoints
      index2=index1+nfocus
      index3=index2+npoints
      top=index3+nfocus
      if(top.gt.canget) then
         write(iout,*) 'need',top,' integer words; have',canget
         call lnkerr('need too much core in m621')
      endif

      call bsj(zc(x),zc(y),zc(z),zc(q),zc(qd),zc(qq),
     $         zc(xcnt),zc(ycnt),zc(zcnt),zc(rad),
     $         zc(xobs),zc(yobs),zc(zobs),epsilon,npoints,nsrcs,
     $         ncnt,nobs,prnt,aben,deltag2,pomf2,self,
     $         nfocus,ac(mark),ac(nark),
     $         zc(phi),zc(ex),zc(ey),zc(ez),zc(efgxx),zc(efgxy),
     $         zc(efgxz),zc(efgyy),zc(efgyz),zc(efgzz),
     $         zc(xq),zc(yq),zc(zq),ac(indx),zc(u),
     $         zc(a2),zc(b2),ncharge,zc(xyzch),zc(charge),
     $         ac(ifocus),ac(icindx),ac(ipindx),ac(index1),ac(index2),
     $         ac(index3),toang,toev,tokcal)
c
c     --- store surface charges and coordinates ---
      call iosys('write integer "number of solvent surface points"'
     $           //' to rwf',1,ncharge,0,' ')
      call iosys('write real "solvent surface coordinates" to rwf',
     $           3*ncharge,zc(xyzch),0,' ')
      call iosys('write real "solvent surface charges" to rwf',
     $           ncharge,zc(charge),0,' ')
c
c     --- store solvent free energy, potential of mean force, and
c         solvent self-energy
      call iosys('write real "solvent free energy" to rwf',
     $           1,deltag2,0,' ')
      call iosys('write real "solvent potential of mean force" to rwf',
     $           1,pomf2,0,' ')
      call iosys('write real "solvent self energy" to rwf',
     $           1,self,0,' ')
c
c     --- if we have converged, override the jump ---
  200 continue
      if(cnvg) then
         call chainx(0)
      else
         call chainx(1)
      endif
c
c
      return
      end
