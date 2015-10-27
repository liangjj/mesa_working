*deck @(#)m202.f	1.3  9/3/91
      program m202
c***begin prologue     m202
c***date written       850601   (yymmdd)
c***revision date      910903   (yymmdd)
c
c   september 3, 1991  bis/rlm at lanl
c            adding a flag which allows us to execute symmetry portions
c            of the code in the presence of ghost or dummy atoms.
c            this is now the default, if sym=(ghosts=off) is present
c            in the options, the code will shut off symmetry in the
c            presence of ghost atoms.
c   7/22/89  bhl at llnl
c            if the option string kohn or no202 is present then
c            this link is not executed
c
c***keywords           m202, link 202, z-matrix, noprtdis, noprtang,
c***                   noangbyz, notetrat, noprtz, nocrowdt, nosymm, coord
c***author             martin, richard (lanl)
c                      binkley, steven, et.al. (gauss82)
c***source             @(#)pm202.f	1.3   9/3/91
c***purpose            determines the coordinates, given the z-matrix,
c                      and analyzes point group symmetry.
c***descriptions
c     this link receives the z-matrix or cartesian coordinates from
c     the rwf files,  determines the point group of the molecule,
c     and produces the standard (3x3) transformation matrices.
c     when input consists of a z-matrix, routine ztoc is called to obtain
c     coordinates.  note that ztoc produces both the regular coordinates
c     and a set in which the dummy atoms have not been deleted.
c
c     m202 currently recognizes the option strings:
c       noprtdis   don't print the distance matrix.
c       noprtang   don't print the interatomic angles.
c       noangbyz   use internal cutoffs instead of the z-matrix to determine
c                  which angles to print.  this is the default, and only
c                  choice, for direct cartesian coordinate input.
c       notetra    don't set angles within 0.001 degree of 109.471
c                  to acos(-1/3).
c       prtz       print the z-matrix and resulting coordinates.
c       nocrowdt   don't abort the run if two atoms are less than
c                  0.5 angstroms apart.
c       nosymm     unconditionally turn symmetry off.
c                  note that symm is still called, and will determine the
c                  framework group, but the molecule is not oriented.
c       noprtsym   turn symmetry printing off.
c       coord      the coordinates were read in cartesian form.
c       timing     print timing statistics for this link.
c
c***references
c
c***routines called
c     m202
c       drum(mdutil)
c       iosys(io)
c       traktm(mdutil)
c       getscm(mdutil)
c       rzmat(local)
c       subvar(util)
c       zprint(util)
c       ztoc(local)
c       wcoord(util)
c       wzmat(local)
c       corpr1(util)
c       dismat(local)
c       angle1(local)
c       angle2(local)
c       rcoord(util)
c       symm(symm)
c       deornt(local)
c       outrep(local)
c       omega(local)
c       rotcon(local)
c       chainx(mdutil)
c
c***end prologue       m202
c
      implicit real*8 (a-h,o-z)
      real*8 z
      integer a
      pointer (p,z(1)), (p,a(1))
      parameter(maxatm=2000,one=1.0d+00)
      logical prtdis,prtang,angbyz,tetrat,prtz,crowdt,cinz,prtrot
      logical prtsym,prtc,inau,logkey,usesym,dmpsym,redsym,norot
      integer iadtwp,wpadti
      integer atmass
      character*4096 ops,info
      character*16 chdum(2*maxatm)
      character*16 chrkey
      character*8 prtflg,symflg,ians
      character pgrp*4,subgrp*4,forgrp*4,forax*1
      real*8 tooclo,toang
      real*8 trvec,rotmat
      common/arcinf/info
      common/io/inp,iout
      dimension trvec(3),rotmat(3,3)
      data prtdis/.true./, prtang/.true./, angbyz/.true./
      data prtc/.true./
      data tetrat/.true./, prtz/.false./, crowdt/.true./
      data cinz/.true./, prtrot/.true./, prtsym/.false./
      data dmpsym/.false./, redsym/.true./, norot/.false./
      data tooclo/0.5d+00/
      data maxop/48/
      data nodum/0/
c
 1000 format(' symmetry turned off by external request.')
 1005 format(' symmetry turned off:')
 1010 format(' atomic calculation.')
 1020 format(' cannot cope with dummy or (shiver) ghost atoms.')
 1030 format(' the molecule is an accidental spherical top.')
 1040 format(' the molecule is an accidental symmetric top.')
 2002 format(1x,'z-matrix orientation:')
 2003 format(1x,'distance matrix (angstroms):')
 2004 format(1x,'input orientation:')
 2005 format(' problem with the z-matrix in link 202.')
 2006 format(' problem with the distance matrix in link 202.')
 2007 format(3d20.12)
 2008 format(i2,3d20.12)
 2010 format(1x,'standard orientation:')
c
      call drum
c      call manmem(0,idum,idum,'m202',0)
c     set up the local options.
      call iosys('read character options from rwf',-1,0,0,ops)
c
c     note...a possible return path
      if(logkey(ops,'no202',.false.,' '))  call chainx(0)
c
      prtdis=logkey(ops,'print=dis',.true.,' ')
      prtang=logkey(ops,'print=ang',.true.,' ')
      prtz=logkey(ops,'print=z',.false.,' ')
      prtc=logkey(ops,'print=c',.false.,' ')
      prtsym=logkey(ops,'print=sym',.false.,' ')
      angbyz=logkey(ops,'geom=angbyz',.true.,' ')
      tetrat=logkey(ops,'geom=tetra',.true.,' ')
      crowdt=logkey(ops,'geom=crowd',.true.,' ')
      usesym=.not.logkey(ops,'sym=off',.false.,' ')
      dmpsym=logkey(ops,'sym=dump',.false.,' ')
      redsym=.not.logkey(ops,'sym=full',.false.,' ')
      forgrp=chrkey(ops,'sym=force_group','    ',' ')
      forax=chrkey(ops,'sym=force_axis',' ',' ')
      norot=logkey(ops,'sym=norotate',.false.,' ')
      cinz=.not.logkey(ops,'geom=coord',.false.,' ')
c
cmp
c
      call iosys('does d2e_coord exist on rwf',0,0,0,ians)
      if(ians.eq.'yes') then
         call iosys('read integer d2e_coord from rwf',1,iicoord,0,' ')
         if (iicoord.eq.1) cinz=.false.
      endif
ccd
      call iosys('does d2e_cycle exist on rwf',0,0,0,ians)
c
c     if print has been turned off elsewhere, e.g. in the course
c     of a geometry optimization, reset the flags.
c
      call iosys('read character "print flag" from rwf',-1,0,0,prtflg)
      if(prtflg.eq.'minimum') then
         prtdis=.false.
         prtang=.false.
         prtrot=.false.
         prtsym=.false.
         prtc=.false.
      endif
c
c     retrieve molecular parameters.
      call iosys('read integer "number of atoms" from rwf',
     $     1,natoms,0,' ')
      call iosys('read integer "number of z-matrix entries" from rwf',
     $     1,nz,0,' ')
      call iosys('read integer "number of variable coordinates" '//
     $     'from rwf',1,nvar,0,' ')
c
c     allocate blank common.
c
c     beginning location        description
c                               integer atomic numbers(natoms).
      kian=1
c                               atomic charges(natoms).
      katchg=iadtwp((kian+max(natoms,nz)))
c                               coordinates(natoms).
      kc=katchg+max(natoms,nz)
c                               atomic masses(natoms).
      atmass=kc+3*max(natoms,nz)
c                               integer z-matrix atomic numbers(nz).
      kianz = wpadti((atmass+natoms))
c                               integer components of the z-matrix(4*nz).
      kiz=kianz+nz
c                               z-matrix bond lengths(nz).
      kbl=iadtwp(kiz+4*nz)
c                               z-matrix bond angles(nz).
      kalpha=kbl+nz
c                               z-matrix dihedral angles(nz).
      kbeta=kalpha+nz
c                               integer z-matrix bond length map(nz).
      klbl=wpadti(kbeta+nz)
c                               integer z-matrix angle map(nz).
      klalph=klbl+nz
c                               integer z-matrix dihedral angle map(nz).
      klbeta=klalph+nz
c                               scratch(8*nz+nz**2)
      iscr1=iadtwp(klbeta+nz)
      iscr2=iscr1+nz
      iscr3=iscr2+nz
      iscr4=iscr3+nz
      iscr5=iscr4+nz
      iscr6=iscr5+nz
      iscr7=iscr6+3*nz
      iscr7=max(iscr7,iscr1+3*natoms)
      iend=wpadti(max(iscr7+nz*nz,iscr7+natoms*natoms))
      maxap3 = natoms + 3
      icnew  = iadtwp(kianz)
      ictmp = icnew + 3*maxap3
      need = wpadti(ictmp + 3*natoms)
      need = max(need,iend)
c      call manmem(need,p,ngot,'geom1',0)
      call getmem(need,p,ngot,'geom1',0)
c
c                                character strings.
c                                names of the z-matrix centers.
      kname=1
c                                basis set names.
      kbtyp=kname+max(nz,natoms)
c
c     skip the z-matrix handling section if direct coordinate input was
c     used.
c
      call iosys('read real angstrom/bohr from rwf',1,toang,0,' ')
      if(cinz) then
c
c        recover the z-matrix variables from the rwf.
         call rzmat('rwf',nz,nvar,a(kianz),a(kiz),z(kbl),
     $              z(kalpha),z(kbeta),a(klbl),a(klalph),
     $              a(klbeta))
c
         if(nvar.ne.0) then
c           recover the substitution values from the rwf.
            call iosys('read real zvalues from rwf',
     $                  -1,z(iscr1),0,' ')
c
c           process possible z-matrix substitutions.
            call subvar(z(kbl),z(kalpha),z(kbeta),a(klbl),
     $                  a(klalph),a(klbeta),z(iscr1),nz,nvar)
         endif
c
c        possibly print the z-matrix, and get the coordinates.
         if(prtz) then
            call iosys('read character z-names from rwf',-1,0,0,
     #                 chdum(kname))
            call iosys('read character "z-matrix basis sets" from rwf',
     $           -1,0,0,chdum(kbtyp))
            call zprint(nz,a(kianz),a(kiz),z(kbl),
     $                  z(kalpha),z(kbeta),a(klbl),
     $                  a(klalph),a(klbeta),chdum(kbtyp),
     $                  chdum(kname))
         endif
         call ztoc(nz,a(kianz),a(kiz),z(kbl),
     $             z(kalpha),z(kbeta),tetrat,numtet,natoms,
     $             a(kian),z(katchg),z(kc),z(iscr6),
     $             z(iscr1),z(iscr2),z(iscr3),z(iscr4),
     $             z(iscr5))
c
c        save full copy of coordinates on rwf.
         call iosys('write real coords&dummies on rwf',
     $              3*nz,z(iscr6),0,' ')
c
c        save z-matrix if tetrahedral angles have been changed,
c        or variables have been substituted.
         if(numtet.ne.0.or.nvar.ne.0) then
            call wzmat('rwf',nz,nvar,a(kianz),a(kiz),z(kbl),
     $                 z(kalpha),z(kbeta),a(klbl),
     $                 a(klalph),a(klbeta))
         endif
c
c        print z-matrix?
         if(prtz) then
            write(iout,2002)
            call corpr1(iout,nz,a(kianz),z(iscr6),toang)
         endif
c
c        compute the distance matrix.
         call dismat(nz,z(iscr6),z(iscr7),toang)
c        test for small distances?
         if(crowdt) then
            call crowd(nz,z(iscr7),a(kianz),tooclo*toang)
         endif
c        print the distance matrix and angles?
         if(natoms.le.1) prtdis=.false.
         if(prtdis) then
            write(iout,2003)
            call iosys('read character z-names from rwf',
     $               -1,0,0,chdum(kname))
            call wlmat(z(iscr7),nz,nz,chdum(kname),
     $                 chdum(kname))
         endif
         if(prtang) then
            if(angbyz) then
               call angle1(iout,nz,a(kiz),a(kianz),z(iscr6))
            else
               call angle2(iout,natoms,nodum,a(kian),z(kc))
            endif
         endif
c
c
      else
c        direct cartesian coordinates input.
c        get the atomic numbers, charges, and coordinates from rwf.
         call rcoord('rwf',natoms,a(kian),z(katchg),z(kc))
	 nodum=0
	 do 10 noat=1,natoms
	    if(z(katchg+noat-1).eq.0.d+0) then
	       nodum=nodum+1
            endif
   10    continue
c
c        print the coordinates?
         if(prtz) then
            write(iout,2004)
            call corpr1(iout,natoms,a(kian),z(kc),toang)
         endif
c
c        compute the distance matrix.
         call dismat(natoms,z(kc),z(iscr7),toang)
c        test for small distances?
         if(crowdt) then
            call crowd(natoms,z(iscr7),a(kian),tooclo*toang)
         endif
c        print the distance matrix and angles?
         if(natoms.le.1) prtdis=.false.
         if(prtdis) then
            call iosys('read character "z-names w/o dummies" from rwf',
     $               -1,0,0,chdum(kname))
            write(iout,2003)
            call wlmat(z(iscr7),natoms,natoms,chdum(kname),
     $                 chdum(kname))
         endif
         if(prtang) call angle2(iout,natoms,nodum,a(kian),z(kc))
      endif
c
c     check to see if symmetry information is requested.
      if (usesym) then
c
c        determine molecular symmetry, check for reorientation if
c        requested, and save symmetry information on the rwf.
c
c        we no longer need the space reserved for the z-matrix.
c        reallocate the distribution of z.
c
c
c        note that at this stage we should always be in atomic units.
         inau=.true.
         call symm(natoms,maxap3,maxop,inau,dmpsym,prtsym,
     $             redsym,forgrp,forax,norot,toang,a(kian),
     $             z(katchg),z(kc),symflg,trvec,rotmat,
     $             pgrp,subgrp,z(icnew),z(ictmp))
         if(symflg.eq.'ok') then
            call deornt(rotmat,1)
c           if(idump.gt.0) call dmprep(natoms,z(iprmut))
cps         call outrep(natoms,z(iprmut))
         else if(symflg.eq.'ghosts') then
            if(logkey(ops,'sym=ghosts=off',.false.,' ')) then
               write(iout,1005)
               usesym=.false.
            else
               call deornt(rotmat,1)
            end if
         else 
            usesym=.false.
            write(iout,1005)
            if(symflg.eq.'atomic') write(iout,1010)
            if(symflg.eq.'axsphtop') write(iout,1030)
            if(symflg.eq.'axsymtop') write(iout,1040)
         end if
      else
         usesym=.false.
cmp         write(iout,1000)
         if(ians.ne.'yes') write(iout,1000)
         call rzero(rotmat,3**2)
         do 20 i=1,3
            rotmat(i,i)=one
   20    continue
      endif
c
c     save the symmetry information.
      call iosys('write character symflg on rwf',len(symflg),
     $            0,0,symflg)
      call iosys('write integer usesym on rwf',1,
     $           usesym,0,' ')
      call iosys('write character "point group" on rwf',len(pgrp),
     $            0,0,pgrp)
      call iosys('write real "rotation matrix" to rwf',9,rotmat,0,' ')
      call omega(natoms,rotmat,trvec,z(kc),z(iscr1),a(kian),
     $           toang,usesym,prtsym,pgrp)
c
c     call omega(maxap3,z(ilmolf),z(imolf),z(ilfwg),z(ifwg),
c    $           isymm,numop,z(itrans),z(inperm),maxop,natoms,
c    $           iprint,idump,trvec,rotmat,z(ngrp),a(kian),
c    $           z(icscr1),toang,z(kc),z(iscr1),z(iscr2),
c    $           z(iscr2),prtsym)
c
c     print the coordinates.
      if(prtflg.ne.'minimum') then
         write(iout,2010)
         call corpr1(iout,natoms,a(kian),z(kc),toang)
      endif
c
c     determine the rotational constants.
      call rotcon(iout,natoms,a(kian),z(kc),pgrp,chdum(1),
     $            prtrot,z(atmass))
c
c     save atomic numbers,charges, and coordinates on rwf.
      call iosys('write integer "number of atoms" to rwf',
     $     1,natoms,0,' ')
c
c      ----- symmetry package uses floated charges, but m102 may have
c            changed them, so re-get proper charges. change this someday.
c
      call iosys('read real "nuclear charges" from rwf',
     $     -1,z(katchg),0,' ')
      call wcoord('rwf',natoms,a(kian),z(katchg),z(kc))
c      call manmem(-ngot,p,idum,'m202',0)
      call getmem(-ngot,p,idum,'m202',0)
c
c
      call chainx(0)
c
c
      stop
      end







