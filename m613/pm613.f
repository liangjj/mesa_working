*deck @(#)pm613.f	5.4 11/28/95
      subroutine pm613(z,a)
c***begin prologue     pm613.f
c***date written       940304  (yymmdd)
c***revision date      11/28/95
c
c***keywords           m613, link 613, density-functional, poisson
c***author             martin, richard (lanl)
c***source             @(#)pm613.f	5.4 11/28/95
c***purpose            evaluates the coulomb potential given the density.
c***description
c
c
c***references         
c                      A.D. Becke, J. Chem. Phys. 88, 2547 (1988)
c                      A.D. Becke and R.M. Dickson,J.Chem.Phys. 89,2993(1988).
c
c***routines called
c***end prologue       pm613.f
c
      implicit none
c     --- input variables -----
c     --- input arrays (unmodified) ---
c     --- input arrays (scratch) ---
      integer a(*)
      real*8 z(*)
c     --- output arrays ---
c     --- output variables ---
c     --- scratch arrays ---
c     --- local variables ---
      integer nbf,nnp,multip,natoms,nae,nbe
      integer nprim,ncont,ntypes,nbtype
      integer ptprim,noprim,ptcont,nocont,start
      integer nocart,nobf,minmom,maxmom,mintyp
      integer nx,ny,nz
      integer cont,ex,top
      integer inp,iout
      integer ndmat
      integer ian,lenxyz,ngrid,mxgrd
      integer c,grid,wts,nradial,nomega
      integer rmax,lmax
      integer d,dtot
      integer maxnbf,mxcont,maxl
      integer canget,left
      integer vwts,rnuc,amu,pwtx,rr,radii
      integer ptrad,radshls,ngb,gblksiz,dftcore,mxgblk,iat,mxgbsiz
      integer iadtwp,wpadti,itobyt
      integer imaxl,charge,size
      integer nugrd,nuwts,nuvwts,newgrd
      integer jmat,values,valuesi
      integer angsiz
      integer intkey,wptoin
      integer maxatm,bigl,minesz
      logical dograd
      data dograd/.false./
      
      parameter (maxnbf=2000,maxatm=2000)
c
      character*8 prtflg
      character*4096 ops
c
      real*8 fpkey
      real*8 zero,two
      real*8 toang
      parameter (zero=0.0d+00,two=2.0d+00)
      real*8 dencut,dmcut,defcut,toosmall
      parameter (defcut=1.0d-16,toosmall=1.0d-50)
c
      logical timeit,debug
      logical prnt,adjust,usesg1,logkey
      parameter (timeit=.true.)
      parameter (debug=.false.)
c
      common /io/     inp,iout
c
      data prnt/.true./
      save prnt
c
 1000 format(1x,'m613:')
 1010 format(5x,'memory available(bytes)',12x,i9)
 1015 format(5x,'maximum block size',17x,i9)
 1016 format(5x,'grid integration information:')
 1017 format(8x,'density-matrix cutoff:           ',1pe8.1)
 1018 format(8x,'density cutoff:                  ',1pe8.1)
 1019 format(8x,'voronoi cells are size adjusted')
 1020 format(5x,'all integrals held in core.')
 1030 format(5x,'# integral triangles in core',i4)
 1040 format(5x,'need(bytes)',18x,i9)
 1050 format(5x,'level-shift used:',f5.2)
 1060 format(8x,'grid size; atom',i3,':',16x,i6,2x,i3,' blocks')
 1070 format(5x,80a1)
 1116 format(5x,'exchange-correlation functional: ',a6,'-',a4)
 1230 format(5x,'total becke charges:',/,
     $      (10x,a8,2x,f11.5))
 1240 format(5x,'total becke charges and spins:',/,
     $      (12x,a8,2f11.6))
c
c     --- get max core available in integers
      call getscm(0,z,canget,'m613:',0)
c
c     ----- recover the options string -----
      call iosys('read character options from rwf',-1,0,0,ops)
      call iosys('read real angstrom/bohr from rwf',1,toang,0,' ')
c
c     ----- check on options and set defaults -----
c     --- quadrature options ---
      rmax=intkey(ops,'scf=radgrid',51,' ')
      lmax=intkey(ops,'scf=lebord',23,' ')
      usesg1=logkey(ops,'scf=sg1',.true.,' ')
      if(usesg1) then
c        for now sg1 always uses at least 51
         rmax=max(rmax,51)
      endif
      dmcut=fpkey(ops,'scf=denmat-cutoff',defcut,' ')
      dencut=fpkey(ops,'scf=density-cutoff',toosmall,' ')
      minesz=intkey(ops,'scf=minesz',100,' ')
      mxgblk=intkey(ops,'scf=maxgblk',5,' ')
      adjust=logkey(ops,'scf=adjustcell',.false.,' ')
c
c     ----- has printing been turned off externally? -----
      call iosys('read character "print flag" from rwf',-1,0,0,prtflg)
      if(prtflg.eq.'minimum') prnt=.false.
      prnt=.true.
      if(prnt) then
         write(iout,1000)
         write(iout,1010) canget*itobyt(1)
         write(iout,1016)
         write(iout,1017) dmcut
         write(iout,1018) dencut
         if(adjust) then
            write(iout,1019) 
         endif
      endif
c
c     ----- get some basic information -----
      call iosys('read integer "number of atoms" from rwf',
     $           1,natoms,0,' ')
      call iosys('read integer "number of basis functions" from rwf',
     $           1,nbf,0,' ')
      nnp=(nbf+1)*nbf/2
      call iosys('read integer "spin multiplicity" from rwf',
     $           1,multip,0,' ')
      call iosys('read integer "number of alpha electrons" from rwf',
     $           1,nae,0,' ')
      call iosys('read integer "number of beta electrons" from rwf',
     $           1,nbe,0,' ')
c
      if (nbf.gt.maxnbf) then
         call lnkerr('character core not long enough')
      end if
c
c     ----- get the number of density matrices ---
      call iosys('read integer "number of hf density matrices" '//
     $           'from rwf',1,ndmat,0,' ')
c
c     --- generate grid points and weights ---
      if (lmax.lt.3 .or. lmax.gt.29)
     $     call lnkerr('m613: invalid lebedev order requested')
      nomega=angsiz(lmax)
      mxgrd=(rmax-1)*nomega
c
c     --- allocate some core
      ian=1
      ngrid=ian+natoms
      ptrad=ngrid+natoms
      radshls=ptrad+rmax*natoms
      ngb=radshls+natoms
      gblksiz=ngb+natoms
      c=iadtwp(gblksiz+mxgblk*natoms)
      grid=c+3*natoms
      wts=grid+3*mxgrd*natoms
      vwts=wts+mxgrd*natoms
      rnuc=vwts+mxgrd*natoms
      amu=rnuc+natoms*natoms
      pwtx=amu+natoms*natoms
      rr=pwtx+natoms
      radii=rr+natoms
      top=wpadti(radii+natoms)
      if (top .gt. canget)
     $     call lnkerr('m613: not enough core for mkgrid')
c     trash everything from rnuc on when mkgrid returns      
      top=wpadti(rnuc)
c
c     --- generate the grid ---
      call iosys('read integer "atomic numbers" from rwf',
     $            -1,a(ian),0,' ')
      call iosys('read real coordinates from rwf',-1,z(c),0,' ')
      call mkgrid(z(c),a(ian),z(grid),z(wts),rmax,lmax,nomega,
     $            nradial,natoms,a(ngrid),mxgrd,z(vwts),z(rnuc),z(amu),
     $            z(pwtx),z(rr),adjust,z(radii),usesg1,
     $            a(radshls),a(ptrad))
c
c     --- redefine mxgrd to be the largest atomic grid generated 
c         this is useful since standard grids may have used less
c         than expected from the default radial and angular orders. 
      newgrd=0
      do 100 iat=1,natoms
         newgrd=max(newgrd,a(ngrid+iat-1))
  100 continue
c     --- pack the points and weights arrays into (newgrd,3,natoms),etc.
c         and redefine pointers
      nugrd=iadtwp(top)
      nuwts=nugrd+3*newgrd*natoms
      nuvwts=nuwts+mxgrd*natoms
      if(newgrd.lt.mxgrd) then
         call pakgrd(z(grid),z(wts),z(vwts),mxgrd,
     $               z(nugrd),z(nuwts),z(nuvwts),newgrd,natoms)
         wts=grid+3*newgrd*natoms
         vwts=wts+newgrd*natoms
         top=wpadti(vwts+newgrd*natoms)
         call vmove(z(grid),z(nugrd),newgrd*3*natoms)
         call vmove(z(wts),z(nuwts),newgrd*natoms)
         call vmove(z(vwts),z(nuvwts),newgrd*natoms)
         mxgrd=newgrd
      endif

c
c     --- retrieve information regarding basis set ---
      call basis(natoms,nbf,nprim,ncont,ntypes,nbtype,lenxyz,
     $           mxcont,maxl,
     $           ptprim,noprim,ptcont,nocont,start,nocart,nobf,
     $           minmom,maxmom,mintyp,nx,ny,nz,cont,ex,
     $           top,z,a)
c     --- allocate core for density matrix and properties
      bigl=max(2,maxl)
      imaxl=top
      d=iadtwp(imaxl+natoms)
      dtot=d+ndmat*nnp
      charge=dtot+nnp
      size=charge+ndmat*natoms
      top=wpadti(size+natoms)
      left=canget-top
c
c     --- we need to figure out how big our grid blocks can be.
c         densty's core requirements are, in units of grid blocks
c
c     real*8 arrays
c     dengrid=1,dengrad=3,ga=1,phi=nbf,grad=3*nbf,tmpgwt=1,scr=1
c     integer arrays
c     nzptrs=1
      dftcore=wptoin(7+4*nbf)+1
c      
c     bfgrd need some temp areas to play with. 
      dftcore=dftcore+1+2*mxcont+3*max((maxl+1),2)
c
c     --- finally we know how much we'll have left for kmatrix.
c         we need dftcore integer word matrices of some length mxgbsiz
      mxgbsiz=left/dftcore
c
c     --- let user override, to avoid using all the memory there is and 
c         choking the machine, but only to DECREASE it, not to increase 
c         it
      mxgbsiz=min(intkey(ops,'scf=mxgbsiz',mxgbsiz,' '),mxgbsiz)
      mxgbsiz=min(mxgrd,mxgbsiz)
c
c     --- now that we know the size of a grid block, let's see how 
c         many grid blocks we need for each atom
c         ngb(iatom) will have number of grid blocks for atom iatom, 
c         gblksiz(block,iatom) gives number of grid points for block in 
c         iatom.
      call gblk(natoms,mxgbsiz,mxgblk,a(ngrid),a(ngb),a(gblksiz))
      if(prnt) then
         do 120 iat=1,natoms
            write(iout,1060) iat,a(ngrid+iat-1),a(ngb+iat-1)
 120     continue 
      endif
c
c     --- set up the pointers ---
      jmat=iadtwp(top)
      values=jmat+nnp*ndmat
      valuesi=wpadti(values)
c
c     --- read the density matrices and let's go
c         first one is the closed shell piece, second one open.
      call iosys('read real "hf density matrix" from rwf',
     #            ndmat*nnp,z(d),0,' ')
c         form the total density
      call rzero(z(dtot),nnp)
      call saxpy(nnp,two,z(d),1,z(dtot),1)
      if(ndmat.eq.2) then
         call vadd(z(dtot),z(dtot),z(d+nnp),nnp)
      endif
      call gofish(z(values),z(dtot),nbf,nnp,z(jmat),ndmat,
     $     ndmat,natoms,mxgrd,
     $     dmcut,dencut,a(ngb),
     $     a(gblksiz),mxgblk,mxgbsiz,a(valuesi),a(ian),
     $     z(c),z(ex),z(cont),a(ptprim),a(noprim),a(nocont),
     $     a(ptcont),mxcont,nprim,ntypes,nbtype,ncont,
     $     a(start),a(nocart),a(nobf),a(maxmom),a(minmom),a(mintyp),
     $     a(nx),a(ny),a(nz),z(grid),z(charge),a(imaxl),
     $     dograd,bigl,minesz,ops)
c
c     --- and exit gracefully ---
      call chainx(0)
c
c
      return
      end
