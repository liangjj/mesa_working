*deck @(#)m202.f	5.1 11/6/94
      program m202
c***begin prologue     m202.f
c***date written       850601   (yymmdd)
c***revision date      11/6/94
c
c   march 17, 1992     rlm at lanl
c            adding an input option which limits the number of atoms
c            which are sent to the symmetry routines.  this is a fix
c            to do copper clusters where there are of the order of 1000
c            point charges surrounding a primary cluster of 10-12 atoms.
c            in this case, the symmetry routines take forever to determine
c            the point group.  the input flag is sym=(patoms=n) where
c            n is the numner of atoms to send to the symmetry routines.
c            for cu2o11 surrounded by point charges, e.g., the route 
c            should contain sym=(patoms=13)
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
c***source             @(#)pm202.f	5.1 11/6/94
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
c***references
c
c***routines called
c
c***end prologue       pm202.f
      implicit none
c     --- input variables -----
c     --- input arrays (unmodified) ---
c     --- input arrays (scratch) ---
      integer a, a1
      real*8 z, z1
c     --- output arrays ---
c     --- output variables ---
c     --- scratch arrays ---
c     --- local variables ---
      integer maxatm
      real*8 one
c
      parameter(maxatm=2000,one=1.0d+00)
c
      integer iadtwp,wpadti
      integer atmass,patoms
      integer inp,iout
      integer maxop,iicoord,natoms,nz,nvar,kian,katchg,kc,kianz
      integer kbl,klbl,iscr,iscr1,iscr2,iscr3,iscr4,iscr5,iscr6,iscr7
      integer kiz,kalpha,kbeta,klalph,klbeta,kbtyp,numtet
      integer top,need,ngot(2),kname,maxap3,icnew,icscr,ictmp,intkey
      integer i, nirrep
      logical prtdis,prtang,angbyz,tetrat,prtz,crowdt,cinz,prtrot
      logical prtsym,prtc,inau,logkey,usesym,dmpsym,redsym,norot
      integer idum
      character*4096 ops,info
      character*16 chdum(2*maxatm)
      character*16 chrkey
      character*8 prtflg,symflg,ians
      character pgrp*4,subgrp*4,forgrp*4,forax*1
      real*8 tooclo,toang
      real*8 trvec(3),rotmat(3,3)
c
      data prtdis/.true./, prtang/.true./, angbyz/.true./
      data prtc/.true./
      data tetrat/.true./, prtz/.false./, crowdt/.true./
      data cinz/.true./, prtrot/.true./
      data dmpsym/.false./, redsym/.true./, norot/.false./
      data tooclo/0.5d+00/
      data maxop/48/
      save prtdis,prtang,angbyz,prtc,tetrat,prtz,crowdt
      save cinz,prtrot,prtsym,dmpsym,redsym,norot,tooclo,maxop
c
      common/arcinf/info
      common/io/inp,iout
      pointer (p,z(1)),(p,a(1))
      pointer (p1,z1(1)),(p1,a1(1))
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
 2009 format(5x,'# primary atoms',i13)
 2010 format(1x,'standard orientation:')
c
      call drum
c     --- set up the local options.
      call iosys('read character options from rwf',-1,0,0,ops)
c
c     --- note...a possible return path
      if(logkey(ops,'no202',.false.,' '))  call chainx(0)
c
      prtdis=logkey(ops,'print=dis',.true.,' ')
      prtang=logkey(ops,'print=ang',.true.,' ')
      prtz=logkey(ops,'print=z',.false.,' ')
      prtc=logkey(ops,'print=c',.false.,' ')
      prtsym=logkey(ops,'print=sym',.true.,' ')
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
c     --- see in what form coordinates exist on rwf
      call iosys('does d2e_coord exist on rwf',0,0,0,ians)
      if(ians.eq.'yes') then
         call iosys('read integer d2e_coord from rwf',1,iicoord,0,' ')
         if (iicoord.eq.1) cinz=.false.
      endif
c
c     --- if print has been turned off elsewhere, e.g. in the course
c         of a geometry optimization, reset the flags.
      call iosys('read character "print flag" from rwf',-1,0,0,prtflg)
      if(prtflg.eq.'minimum') then
         prtdis=.false.
         prtang=.false.
         prtrot=.false.
         prtsym=.false.
         prtc=.false.
      endif
c
c     --- retrieve molecular parameters.
      call iosys('read integer "number of atoms" from rwf',
     $            1,natoms,0,' ')
      call iosys('read integer "number of z-matrix entries" from rwf',
     $            1,nz,0,' ')
      call iosys('read integer "number of variable coordinates" '//
     $            'from rwf',1,nvar,0,' ')
c
c     --- allocate memory for arrays used in all of m202.
c     beginning location        description
c                               integer atomic numbers(natoms).
      kian=1
c                               atomic charges(natoms).
      katchg=iadtwp((kian+max(natoms,nz)))
c                               coordinates(natoms).
      kc=katchg+max(natoms,nz)
c                               atomic masses(natoms).
      atmass=kc+3*max(natoms,nz)
c                               scratch
      iscr=atmass+natoms
      need=wpadti(iscr+natoms*natoms)
      call getmem(need,p,ngot(1),'m202',0)
      call iosys('read real angstrom/bohr from rwf',1,toang,0,' ')
c
c     what is done next depends on whether we are doing z matrix or cartesian geometry.
c
      if(cinz) then
c
c        ok, we want a z matrix
c
c                                character strings.
c                                names of the z-matrix centers.
         kname=1
c                                basis set names.
         kbtyp=kname+max(nz,natoms)
c
c
c        allocate the core for the z matrix subroutines.
c
c                               integer z-matrix atomic numbers(nz).
         kianz=1
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
         top=wpadti(max(iscr7+nz*nz,iscr7+natoms*natoms))
c        call getscm(top,a(kian),ngot,'geom1',0)
         call getmem(need,p1,ngot(2),'z-mat',0)

c
c        ---  recover the z-matrix variables from the rwf.
         call rzmat('rwf',nz,nvar,a1(kianz),a1(kiz),z1(kbl),
     $              z1(kalpha),z1(kbeta),a1(klbl),a1(klalph),
     $              a1(klbeta))
c
         if(nvar.ne.0) then
c           --- recover the substitution values from the rwf.
            call iosys('read real zvalues from rwf',
     $                  -1,z1(iscr1),0,' ')
c
c           --- process possible z-matrix substitutions.
            call subvar(z1(kbl),z1(kalpha),z1(kbeta),a1(klbl),
     $                  a1(klalph),a1(klbeta),z1(iscr1),nz,nvar)
         endif
c
c        --- possibly print the z-matrix, and get the coordinates.
         if(prtz) then
            call iosys('read character z-names from rwf',-1,0,0,
     $                 chdum(kname))
            call iosys('read character "z-matrix basis sets" from rwf',
     $                 -1,0,0,chdum(kbtyp))
            call zprint(nz,a1(kianz),a1(kiz),z1(kbl),
     $                  z1(kalpha),z1(kbeta),a1(klbl),
     $                  a1(klalph),a1(klbeta),chdum(kbtyp),
     $                  chdum(kname))
         endif
         call ztoc(nz,a1(kianz),a1(kiz),z1(kbl),
     $             z1(kalpha),z1(kbeta),tetrat,numtet,natoms,
     $             a(kian),z(katchg),z(kc),z1(iscr6),
     $             z1(iscr1),z1(iscr2),z1(iscr3),z1(iscr4),
     $             z1(iscr5))
c
c        --- save full copy of coordinates on rwf.
         call iosys('write real coords&dummies on rwf',
     $              3*nz,z1(iscr6),0,' ')
c
c        --- save z-matrix if tetrahedral angles have been changed,
c            or variables have been substituted.
         if(numtet.ne.0.or.nvar.ne.0) then
            call wzmat('rwf',nz,nvar,a1(kianz),a1(kiz),z1(kbl),
     $                 z1(kalpha),z1(kbeta),a1(klbl),
     $                 a1(klalph),a1(klbeta))
         endif
c
c        --- print z-matrix?
         if(prtz) then
            write(iout,2002)
            call corpr1(iout,nz,a1(kianz),z1(iscr6),toang)
         endif
c
c        --- compute the distance matrix.
         call dismat(nz,z1(iscr6),z1(iscr7),toang)
c        --- test for small distances?
         if(crowdt) then
            call crowd(nz,z1(iscr7),a1(kianz),tooclo*toang)
         endif
c        --- print the distance matrix and angles?
         if(natoms.le.1) prtdis=.false.
         if(prtdis) then
            write(iout,2003)
            call iosys('read character z-names from rwf',
     $               -1,0,0,chdum(kname))
            call wlmat(z1(iscr7),nz,nz,chdum(kname),
     $                 chdum(kname))
         endif
         if(prtang) then
            if(angbyz) then
               call angle1(iout,nz,a1(kiz),a1(kianz),z1(iscr6))
            else
               call angle2(iout,natoms,a(kian),z(kc))
            endif
         endif
c      
c        recover memory from z-matrix step
c
         call getmem(-ngot(2),p1,idum,'z-mat',idum)
c
c
      else
c        --- direct cartesian coordinates input.
c            get the atomic numbers, charges, and coordinates from rwf.
         call rcoord('rwf',natoms,a(kian),z(katchg),z(kc))
c
c        --- print the coordinates?
         if(prtz) then
            write(iout,2004)
            call corpr1(iout,natoms,a(kian),z(kc),toang)
         endif
c
c        --- compute the distance matrix.
         call dismat(natoms,z(kc),z(iscr),toang)
c        --- test for small distances?
         if(crowdt) then
            call crowd(natoms,z(iscr),a(kian),tooclo*toang)
         endif
c        --- print the distance matrix and angles?
         if(natoms.le.1) prtdis=.false.
         if(prtdis) then
            kname=1
            call iosys('read character "z-names w/o dummies" from rwf',
     $               -1,0,0,chdum(kname))
            write(iout,2003)
            call wlmat(z(iscr),natoms,natoms,chdum(kname),
     $                 chdum(kname))
         endif
         if(prtang) call angle2(iout,natoms,a(kian),z(kc))
      endif
c
c     --- check to see if symmetry information is requested.
      pgrp='none'
      if (usesym) then
c
c        ---determine molecular symmetry, check for reorientation if
c           requested, and save symmetry information on the rwf.
c
c           we no longer need the space reserved for the z-matrix.
c           reallocate the distribution of z.
c
         maxap3 = natoms + 3
         icnew  = 1
         ictmp = icnew + 3*maxap3
         top = wpadti(ictmp + 3*natoms)
c         call getscm(top,z(icnew),ngot,'geom2',0)
         call getmem(top,p1,ngot(2),'sym',0)
c
c        --- note that at this stage we should always be in atomic units.
         inau=.true.
         patoms=intkey(ops,'sym=patoms',natoms,' ')
         if(patoms.ne.natoms) then
            write(iout,2009) patoms
         endif
         call symm(patoms,maxap3,maxop,inau,dmpsym,prtsym,
     $             redsym,forgrp,forax,norot,toang,a(kian),
     $             z(katchg),z(kc),symflg,trvec,rotmat,
     $             pgrp,subgrp,z1(icnew),z1(ictmp))
         if(symflg.eq.'ok') then
            call deornt(rotmat,1)
c           if(idump.gt.0) call dmprep(natoms,z(iprmut))
cps         call outrep(natoms,z(iprmut))
         else if(symflg.eq.'atomic') then
            call deornt(rotmat,1)
            write(iout,1010)
         else if(symflg.eq.'ghosts') then
            if(logkey(ops,'sym=ghosts=off',.false.,' ')) then
               write(iout,1005)
               usesym=.false.
            else
               call deornt(rotmat,1)
c               write(iout,1020)
            end if
         else
            usesym=.false.
            write(iout,1005)
            if(symflg.eq.'axsphtop') write(iout,1030)
            if(symflg.eq.'axsymtop') write(iout,1040)
         end if
         call getmem(-ngot(2),p1,idum,'sym',idum)
      else
         usesym=.false.
         if(prtsym) then
            write(iout,1000)
         endif
         call rzero(rotmat,3**2)
         do 20 i=1,3
            rotmat(i,i)=one
   20    continue
      endif
c
c     --- save the symmetry information.
      call iosys('write character symflg on rwf',len(symflg),
     $            0,0,symflg)
      call iosys('write integer usesym on rwf',1,
     $           usesym,0,' ')
      call iosys('write character "point group" on rwf',len(pgrp),
     $            0,0,pgrp)
      call iosys('write real "rotation matrix" to rwf',9,rotmat,0,' ')
      call omega(natoms,rotmat,trvec,z(kc),z(iscr),a(kian),
     $           toang,usesym,dmpsym,pgrp)
c
c     --- print the coordinates.
      if(prtflg.ne.'minimum') then
         write(iout,2010)
         call corpr1(iout,natoms,a(kian),z(kc),toang)
      endif
c
c     --- determine the rotational constants.
      call rotcon(iout,natoms,a(kian),z(kc),pgrp,chdum(1),
     $            prtrot,z(atmass))
c
c     --- save atomic numbers,charges, and coordinates on rwf.
      call iosys('write integer "number of atoms" to rwf',
     $           1,natoms,0,' ')
c
c     --- symmetry package uses floated charges, but m102 may have
c         changed them, so re-get proper charges. change this someday.
      call iosys('read real "nuclear charges" from rwf',
     $           -1,z(katchg),0,' ')
      call wcoord('rwf',natoms,a(kian),z(katchg),z(kc))
      call getmem(-ngot(1),p,ngot(1),'m202',0)
c
c
      call chainx(0)
c
c
      stop
      end
