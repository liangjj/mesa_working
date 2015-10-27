*deck @(#)pm101.f	5.2  2/5/95
      program m101
c***begin prologue     pm101.f
c***date written       850601   (yymmdd)
c***revision date      931108   (yymmdd)
c
c    23 january 1995   rlm at lanl
c       adding capability to read grid in the $geom section
c    27 october 1994   rlm at lanl
c       adding capability to read solvent radii in the $geom section
c
c    8  november 1993 rlm at lanl
c       adding capability to read an initial guess at the full second
c       derivative matrix from the input stream.
c
c    14 december 1988 bhl at llnl
c       unicos changes for memory allocation
c
c    11 january 1988
c       page fix for geom=(rdchk)
c
c     2 september 1987  pws at lanl
c        changing the core alloaction to use integer and real arrays,
c        as this will hopefully fix some problems associated with
c        dummy atoms.
c
c     1 december 1986   pws at lanl
c        changing 'namchk' and iosys open to character variables
c
c        modified 22 july 1986 at lanl by pws
c                 writing out the atom names [chdum(kname)] as they are
c                 input to the z-matrix.
c
c***keywords
c                      m101, link 101, input,
c                      geometry, basis set, $geom,
c                      coord, rdchk, inau, inrad
c***author             martin, richard  (lanl)
c                      binkley, steven  (gauss82)
c***source             @(#)pm101.f	5.2   2/5/95
c***purpose
c                      reads the geometry, and basis set information.
c***description
c     m101 reads the geometry specifications, either in the
c     form of a symbolic z-matrix with associated variables and constants,
c     or a direct cartesian coordinate list.
c     for additional information see the documentation for $geom
c
c     options:
c       coord     the coordinates are to be read in cartesian form.
c       rdchk     read the geometry from the checkpoint file.
c       inau      the distances are in atomic units.
c       inrad     the angles are in radians.
c       print     print geometry information.
c
c       input:    data on cards in the input deck and options in /ops/.
c       output:   in atomic units (bohrs/radians).
c                 nz     ... the number of cards in the z-matrix.
c                 nvar   ... the number of variables in the z-matrix,
c                            equal to the number of symbols defined
c                            variables section of the input.
c                 the z-matrix information is contained in several arrays.
c                 ianz   ... atomic numbers; integer(nz).
c                 iz     ... connectivity data; integer(4,nz).
c                 bl     ... bond lengths; real(nz).
c                 alpha  ... valence angles; real(nz).
c                 beta   ... dihedral or second valence angle; real(nz).
c                 lbl    ... used to map from the array of
c                            variables in the substitution information
c                            to the array bl; integer(nz).
c                            if lbl(n)=0, no substitution is required
c                               to get bl(n), it already contains the
c                               proper value.
c                            if lbl(n)= i, bl(n)=values(i).
c                            if lbl(n)=-i, bl(n)=-values(i).
c                 lalpha ... analogous to lbl for alpha; integer(nz).
c                 lbeta  ... analogous to lbl for beta; integer(nz).
c
c                 z-matrix substitutions are contained in the arrays:
c                 anames ... alphanumeric names of the variables;
c                            character*16 (nvar).
c                 values ... the corresponding numeric values. these
c                            are altered in the course of a geometry
c                            optimization; real(nvar).
c                 intvec ... an array of the integer values following
c                            the symbol and value on a line of the
c                            variables section; integer(nvar).
c                 fpvec  ... a corresponding array containing the floating
c                            point number. the use of these two vectors
c                            depend on the route; real(nvar).
c
c***references
c
c***routines called
c
c***end prologue       pm101.f
      implicit none
c     --- input variables -----
c     --- input arrays (unmodified) ---
c     --- input arrays (scratch) ---
      integer a
      real*8 z
c     --- output arrays ---
c     --- output variables ---
c     --- scratch arrays ---
c     --- local variables ---
      integer maxatm,maxchr
c
      parameter (maxatm=2000)
      parameter (maxchr=12*maxatm)
c
      integer iadtwp,wpadti
      integer inp,iout
      integer nel
      integer ian,atchg,c,ianz,iz,bl,alpha,beta,lbl,lalpha,lbeta,value
      integer intv,fpvec,scr1,scr2,scr3,scr4,scr5,scr6
      integer need,btyp,bname,gtyp,gname,sym,cntr,ctop
      integer nz,natoms,nvar,iat,i,j
      integer name,tetrat,numtet,ftric,nvvc,top
      integer radius
      integer ngot, idum
      logical logkey,cinz,rdinp,inau,inrad,positn
      character*16 chdum(maxchr)
      character*16 defbas,defgrd,chrkey,srcd2e
      character*4096 ops,info
      character*2 el(106)
      character*128 namchk
      character card*80
      real*8 toang
      pointer(p,z(1)),(p,a(1))
c
      data nel/92/
      data cinz/.true./, rdinp/.true./ 
      data inau/.false./, inrad/.false./
      save nel,cinz,rdinp,inau,inrad
c
      common/arcinf/info
      common/io/inp,iout
c
 1000 format(/,1x,'m101:  process geometry input.')
 1050 format('     generated z-matrix:')
 1070 format('     cartesian coordinates read from the input stream:')
 1080 format(8f10.6)
c
      call drum
c     --- announce our presence
      write(iout,1000)
c
c     --- set up local options.
      call iosys('read character options from rwf',0,0,0,ops)
      if(logkey(ops,'geom=coord',.false.,' ')) cinz=.false.
      if(logkey(ops,'geom=rdchk',.false.,' ')) rdinp=.false.
      if(logkey(ops,'geom=chk',.false.,' ')) rdinp=.false.
      if(logkey(ops,'geom=inau',.false.,' ')) inau=.true.
      if(logkey(ops,'geom=inrad',.false.,' ')) inrad=.true.
      srcd2e=chrkey(ops,'opt=srcd2e',' ',' ')
c
c     --- allocate blank common.
c     beginning location   description
c                          integer atomic numbers(natoms).
      ian=1
c                          atomic charges(natoms).
      atchg=iadtwp(ian+maxatm)
c                          coordinates(3*maxatm).
      c=atchg+maxatm
c                          solvent radii(maxatm)
      radius=c+3*maxatm
c
c                          z-matrix section.
c                          z-matrix atomic numbers(maxatm).
      ianz=wpadti(radius+maxatm)
c                          integer components of the z-matrix(4*maxatm).
      iz=ianz+maxatm
c                          z-matrix bond lengths(maxatm).
      bl=iadtwp(iz+4*maxatm)
c                          z-matrix bond angles(maxatm).
      alpha=bl+maxatm
c                          z-matrix dihedral angles(maxatm).
      beta=alpha+maxatm
c                          z-matrix bond length map(maxatm).
      lbl=wpadti(beta+maxatm)
c                          z-matrix angle map(maxatm).
      lalpha=lbl+maxatm
c                          z-matrix dihedral angle map(maxatm).
      lbeta=lalpha+maxatm
c                            variables section.
c                            values of the variables.
      value=iadtwp(lbeta+maxatm)
c                            integer values.
      intv=wpadti(value+maxatm)
c                            floating point values.
      fpvec=iadtwp(intv+maxatm)

      scr1=fpvec+maxatm
      scr2=scr1+3*maxatm
      scr3=scr2+3*maxatm
      scr4=scr3+3*maxatm
      scr5=scr4+3*maxatm
      scr6=scr5+3*maxatm
      need=wpadti(scr6+3*3*maxatm)
c                            character strings.
c                            names of the variables.
      name=1
c                             basis set and grid names.
      btyp=name+maxatm
      bname=btyp+maxatm
c
      gtyp=bname+maxatm
      gname=gtyp+maxatm
c                             scratch.
      sym=gname+maxatm
      cntr=sym+5*maxatm
c
      ctop=cntr+maxatm-1
      if (ctop.gt.maxchr) then
         call lnkerr('not enough character core in m101')
      end if
c
c     --- don't need the z-matrix if we're using cartesians instead.
      if(.not.cinz) need=ianz
c
c      call getscm(need,a,maxcor,'m101',0)
      call getmem(need,p,ngot,'m101',0)
c
c
      call fillel(0,nel,el)
      if(rdinp) then
         if(.not.positn('$geom',card,inp))
     $      call lnkerr('no $geom section found.')
      endif
      nz = 0
      natoms = 0
      nvar = 0
      call iosys('read real angstrom/bohr from rwf',1,toang,0,' ')
c
c     --- let's go.
      if(rdinp) then
         if(cinz) then
c           read a symbolic z-matrix from input stream.
            call zget(nz,nvar,maxatm,toang,inau,inrad,a(ianz),
     $           z(atchg),a(iz),z(bl),z(alpha),
     $           z(beta),a(lbl),a(lalpha),a(lbeta),
     $           z(value),z(fpvec),a(intv),
     $           chdum(name),chdum(sym),chdum(cntr),chdum(btyp),
     $           chdum(gtyp),z(radius))
c
         else
c           read cartesian coordinates from input stream.
            call coord(natoms,a(ian),z(atchg),z(c),z(radius),
     $                 chdum(bname),chdum(gname),chdum(cntr),el,inau,
     $                 toang)
         endif
      else
c        --- open the checkpoint file.
         call iosys('read character "checkpoint filename" from rwf',
     $               0,0,0,namchk)
         call iosys('open chk as old',1000000,0,0,namchk)
         call iosys('read integer "number of atoms" from chk',1,
     $               natoms,0,' ')
         call iosys('read integer "number of z-matrix entries" '//
     $              'from chk',1,nz,0,' ')
         call iosys('read integer "number of variable coordinates" '//
     $              'from chk',1,nvar,0,' ')
         call iosys('read character "atomic basis name" from chk',
     $              -1,0,0,chdum(bname))
         call iosys('read character "atomic grid name" from chk',
     $              -1,0,0,chdum(gname))
         call iosys('read character "z-names w/o dummies" from chk',
     $              -1,0,0,chdum(cntr))
         if(cinz) then
c           read z-matrix from the checkpoint file.
            call rzmat('chk',nz,nvar,a(ianz),a(iz),z(bl),
     $                 z(alpha),z(beta),a(lbl),
     $                 a(lalpha),a(lbeta))
c
            call iosys('read character "z-matrix basis sets" from chk',
     $                 -1,0,0,chdum(btyp))
            call iosys('read character "z-matrix atomic grids"'
     $                //' from chk',-1,0,0,chdum(gtyp))
            do 100 j= 1,nz
               z(atchg-1+j)=float(a(ianz-1+j))
 100        continue
cmp
c
            if(nvar.gt.0) then
c              --- read the substitution data from the check file.
               call rzsub('chk',nvar,chdum(name),z(value),
     $                    a(intv),z(fpvec))
            endif
         else
c           --- read cartesian coordinates from checkpoint file.
            call iosys('read integer "atomic numbers" from chk',
     $                 -1,a(ian),0,' ')
            call iosys('read real coordinates from chk',
     $                 -1,z(c),0,' ')
         endif
c        --- close the checkpoint file.
         call iosys('close chk',namchk,0,0,' ')
      endif
c
c     --- check for nonsense.
      if(nz.le.0.and.natoms.le.0) call lnkerr('there are no atoms!')
c
c     --- store the z-matrix on the rwf and prepare the atomic numbers
c         and charges.
      if(cinz) then
c        --- write the /zmat/.
c
         call wzmat('rwf',nz,nvar,a(ianz),a(iz),z(bl),
     $              z(alpha),z(beta),a(lbl),a(lalpha),
     $              a(lbeta))
c
c
c        --- write the z-matrix names
         call iosys('write character z-names to rwf',16*nz,0,0,
     $               chdum(cntr))
c        --- write /zsubst/
         if(nvar.ne.0) then
            call wzsub('rwf',nvar,chdum(name),z(value),a(intv),
     $                 z(fpvec))
         endif
c
c        --- if not specified in the options,
c            default the basis to double-zeta and the grid to sg1.
         defbas=chrkey(ops,'basis','dz',' ')
         defgrd=chrkey(ops,'grid','sg1',' ')
         do 20 iat=1,nz
            if (a(ianz+iat-1).ge.0) then
               if(chdum(btyp+iat-1).eq.' ') chdum(btyp+iat-1)=defbas
               if(chdum(gtyp+iat-1).eq.' ') chdum(gtyp+iat-1)=defgrd
            end if
 20      continue
c
c        --- write the basis sets associated with z-matrix entries.
         call iosys('write character "z-matrix basis sets" to rwf',
     $        16*nz,0,0,chdum(btyp))
c        --- write the grids sets associated with z-matrix entries.
         call iosys('write character "z-matrix atomic grids" to rwf',
     $        16*nz,0,0,chdum(gtyp))
c
c        --- prepare the atomic numbers and charges.
c            one must distinguish between ghosts and dummies here.
c            dummies have an atomic number equal to -1, and are used to
c            mark a coordinate only.  they are omitted from the list of
c            atoms. ghosts have an atomic number of 0, and may have an 
c            associated basis set and grid.
         natoms=0
         do 10 i=1,nz
c           --- omit dummies.
            if(a(ianz+i-1).ge.0) then
               natoms=natoms+1
               a(ian+natoms-1)=a(ianz+i-1)
               chdum(bname+natoms-1)=chdum(btyp+i-1)
               chdum(gname+natoms-1)=chdum(gtyp+i-1)
               chdum(cntr+natoms-1)=chdum(cntr+i-1)
               z(atchg+natoms-1)=z(atchg+i-1)
               z(radius+natoms-1)=z(radius+i-1)
            endif
  10    continue
c
c       --- generate cartesian coordinates from the z-matrix
        if(nvar.ne.0) then
c          --- process possible z-matrix substitutions.
           call subvar(z(bl),z(alpha),z(beta),a(lbl),
     $                 a(lalpha),a(lbeta),z(value),nz,nvar)
        endif
        call ztoc(nz,a(ianz),a(iz),z(bl),
     $            z(alpha),z(beta),tetrat,numtet,natoms,
     $            a(ian),z(atchg),z(c),z(scr6),
     $            z(scr1),z(scr2),z(scr3),z(scr4),
     $            z(scr5))
c        --- save full copy of coordinates on rwf
         call iosys('write real coords&dummies to rwf',3*nz,
     $               z(scr6),0,' ')
      else
c
c        --- tidy up any atoms using the default basis.
         defbas=chrkey(ops,'basis','dz',' ')
         defgrd=chrkey(ops,'grid','sg1',' ')
         do 30 iat=1,natoms
            if(chdum(bname+iat-1).eq.' ') chdum(bname+iat-1)=defbas
            if(chdum(gname+iat-1).eq.' ') chdum(gname+iat-1)=defgrd
 30      continue
      endif
c
c     --- save often used information in global common block.
      call iosys('write integer "number of atoms" to rwf',
     $            1,natoms,0,' ')
      call iosys('write integer "number of z-matrix entries" to rwf',
     $            1,nz,0,' ')
      call iosys('write integer "number of variable coordinates" '//
     $           'to rwf',1,nvar,0,' ')
c
c     --- write the coordinates ,atomic numbers, basis sets, and
c         solvent radii to the rwf.
      call wcoord('rwf',natoms,a(ian),z(atchg),z(c))
      call iosys('write character "atomic basis name" to rwf',
     $           16*natoms,0,0,chdum(bname))
      call iosys('write character "atomic grid name" to rwf',
     $           16*natoms,0,0,chdum(gname))
      call iosys('write character "z-names w/o dummies" to rwf',
     $           16*natoms,0,0,chdum(cntr))
      call iosys('write real "solvent radii" to rwf',
     $           natoms,z(radius),0,' ')
c
c     --- print the geometry information.
      if(cinz) then
         if(.not.rdinp) then
            write(iout,*) '    z-matrix read from checkpoint file '
         endif
         call subvar(z(bl),z(alpha),z(beta),
     $               a(lbl),a(lalpha),a(lbeta),
     $               z(value),nz,nvar)
c
         call zprint(nz,a(ianz),a(iz),z(bl),
     $               z(alpha),z(beta),a(lbl),
     $               a(lalpha),a(lbeta),chdum(btyp),
     $               chdum(name))
      else
         if(.not.rdinp) then
            write(iout,*) '    cartesian coordinates read from'
     $                  //'checkpoint file '
         endif
         call cprint(natoms,a(ian),z(atchg),z(c),toang,
     $        chdum(bname))
      endif
      call getmem(-ngot,p,idum,'m101',idum)
c
c     --- perhaps read an initial guess at the full second derivative matrix
c         in cartesian coordinates. only the lower triangle is read.
c         note that m732 will convert this into internal coordinates.
c         there is a gradient dependent part to this transformation,however,
c         and so we must also provide m732 with cartesian gradients.
c         trick it here by loading zeros into the gradient.
      if (srcd2e.eq.'rdinp') then
         if(.not.positn('$cfc',card,inp))
     $      call lnkerr('no $cfc section found.')
c        --- make sure we have enough core ---
         ftric=1
         nvvc=3*natoms*(3*natoms+1)/2
         need=wpadti(ftric+nvvc)
         call getmem(need,p,ngot,'m101',0)
c        read(inp,1080) (z(i),i=1,nvvc)
         read(inp,*) (z(i),i=1,nvvc)
         call iosys('write real "cartesian second derivatives" to rwf',
     $              nvvc,z(ftric),0,' ')
c        assume zero gradients for now.
         call rzero(z(ftric),3*natoms)
         call iosys('write real "cartesian first derivatives" to rwf',
     $              3*natoms,z(ftric),0,' ')
         call getmem(-ngot,p,idum,'m101',idum)
      endif
c
c
      call chainx(0)
c
c
      stop
      end
