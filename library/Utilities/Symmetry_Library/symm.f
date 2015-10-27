*deck @(#)symm.f	5.1  11/6/94
      subroutine symm(natoms,maxap3,maxop,inau,dump,prtsym,
     $                redsym,forgrp,forax,norot,toang,ian,atmchg,cin,
     $                symflg,trvec,rotmat,pgrp,subgrp,cnew,ctmp)
c***begin prologue     symm
c***date written       850601  (yymmdd)
c***revision date      930630  (yymmdd)
c      june 30,  1993  rlm at lanl
c                        adding capability to handle atoms in d2h.
c      march 23, 1993  rlm at lanl
c                      when sym=norotate is specified, set translation
c                      vector to zero.
c      march 17,1990   rlm at lanl
c                      extensive modifications to oper and the routines
c                      which it calls.
c***keywords           symmetry, point group
c***author             defrees, doug (cmu)
c                      saxe,paul and martin,richard (lanl)
c***purpose            main driver for the symmetry package.
c***description
c
c     this is the main driver for the symmetry package.
c
c     given the coordinates and atomic numbers (or any other identifying
c     feature such as atomic weights) this package determines:
c     the molecular point group.
c     the standard orientation of the molecule in cartesian space.
c     the atom permutation array under each operation of the group.
c     the 3x3 rotation matrices and the basis function transformation
c     matrices for the operations of the group.
c     the symmetry adapted linear combinations (salc's) of the basis set
c     which span the irreducible representations.
c
c
c     input arguments:
c     natoms ... number of atoms.
c     maxap3 ... natoms+3 for dimensioning.
c     maxop  ... maximum number of operations.  normally 48; i.e. oh.
c     inau   .. .true.  if coordinates are in bohr.
c               .false. if coordinates are in angstroms.
c     dump   ... logical debugging dump flag.
c     prtsym ... logical print flag.
c     redsym ... logical symmetry reduction flag.
c                .true. (default) work in d2h and subgroups.
c                .false.          work in full symmetry group.
c     forax  ... force the principal rotation axis.  the default is the z-axis.
c                this option is useful in those situations where the standard
c                geometry might not be precisely the alignment you would like.
c                see the discussion under forgrp.
c     forgrp ... force the point group.
c                if non-blank, the user requests the calculation to be
c                specifically done for a specific group given by the character string
c                forgrp. this is usually used when the investigator wants to
c                perform the calculation in a subgroup of the full point group,
c                but does not want to accept the default principal axis definition
c                as the z-axis. for example, suppose the full point group is
c                d3h, with the z-axis the c3 axis. a calculation in c2v could
c                be performed in which the principal axis is defined as c2(y).
c                by specifying norot=.true., forgrp='c2v', and forax='y'.
c     norot  ... rotation to standard orientation?
c                .true.   don't rotate the molecule before symmetry analysis.
c                in this case, the user is responsible for defining
c                coordinate axes compatible with the conventions
c                used for the group generators and basis function
c                transformations. the order and sense of the generators
c                is:   
c                      c2(z)      if it exists.
c                      c2'(y)     if it exists.
c                      sigma(xy)  
c                      inversion  if they exist.
c
c                if the user wishes to override these defaults, he or she
c                can do so by specifying the point group of interest in forgrp,
c                and the principal axis in forax. this override will work only
c                for d2h and subgroups.
c
c                .false.  rotate to standard orientation before symmetry
c                         analysis.
c     toang  ... conversion factor from bohr to angstroms.
c     ian    ... atomic numbers, used for producing atomic symbols.
c     atmchg ... atomic numbers or atomic masses, as appropriate.
c     cin    ... input coordinates dimensioned (3,natoms).
c
c
c     implicit input:
c     common/io/inp,iout   input/output unit numbers.
c     basis set information (see m102)
c     this information is read from the rwf.
c     maxmom ... the maximum angular momentum recognized by the program.
c     nbf    ... the number of basis functions.
c     ncart  ... the number of cartesian components in each angular
c                momentum. ncart(0:maxmom).
c     nx     ... the power of x in the polynomial portion of an
c                angular momentum block.
c     ny     ... the power of y in the polynomial.
c     nz     ... the power of z in the polynomial.
c     momatm ... the maximum momentum of any shell on a particular atom.
c                momatm(natoms).
c     nocont ... the number of contracted functions of each angular
c                momentum on each atom. nocont(natoms,0:maxmom)
c     bfstrt ... the starting index(-1) of each angular momentum group
c                on each atom.
c
c
c     output arguments:
c     symflg ... flag for whether symmetry was successfully determined:
c     'ok'       ... symmetry has successfully been set up.
c     'atomic'   ... atomic calculation.
c     'ghosts'   ... ghost atoms present.
c     'axsphtop' ... the molecule is an accidental spherical top.
c     'axsymtop' ... the molecule is an accidental symmetric top.
c     trvec  ... translation vector for conversion to standard
c                orientation.
c     rotmat ... rotation matrix for conversion to standard orientation.
c     pgrp   ... name of the group, up to 4 characters.
c     subgrp ... name of the subgroup, up to 4 characters.
c     cnew   ... the coordinates for the standard orientation
c                and the rotation matrix which produces them.
c                (maxap3,3).
c
c     implicit output:
c     this information is written to the rwf; see routine oper.
c
c     ngen   ... number of generators.
c     nirrep ... number of irreducible representations.
c     nop    ... number of symmetry operations.
c     trans  ... coordinate transformations: (3,3,nop)
c                trans(,,i) gives the rotation applied to
c                the coordinate system by operation i.
c     bftran ... basis function transformations: (ncart,ncart,nop)
c                bftran(,,i) gives the transformation matrix
c                describing how basis functions mix under
c                operation i. there are separate matrices for
c                each angular momentum block.
c     ptbftr ... pointer array to address bftran: (0:maxmom).
c     atprmt ... atom permutation map:(natoms,nop).
c                atprmt(i,j) is the atom to which atom i
c                is converted by operation j.  
c     symat  ... atom sets:(natoms).
c                which symmetry distinct set does atom i belong to? 
c     ns     ... number of symmetry distinct atoms.
c                this is the number of sets of atoms which transform
c                into one another under the operations of the group.
c     mcu    ... maximum number of symmetry related atoms.
c     nsymat ... number of symmetry related atoms:(ns).
c                the number of atoms in each symmetry related
c                set.
c     relatm ... symmetry related atoms:(mcu*ns).
c                the list of atom indices comprising each set.
c     ptsc   ... symmetry contraction pointers:(0:maxmom,ns)
c                pointer array into the symmetry contraction matrices.
c     sc     ... symmetry contraction matrices.
c     numso  ... number of salc's in each irreducible representation:(nirrep)
c     lambda ... degeneracy of each irreducible representation:(nirrep)
c     lirrep ... labels of irreducible representations:(nirrep)
c     coeffs ... the salc transformation matrix:(nbf,nsalc).
c                nsalc is equal to nbf.
c     naords ... number of ao reduction sets.
c     maxsao ... maximum number of related ao's.
c     aords  ... ao reduction sets:(maxsao,naords)
c
c
c     scratch arrays:
c     ctmp   ... (natoms,3).
c     a      ... beginning of integer blank common.
c     z      ... beginning of real blank common. 
c            ... oper assumes this is the beginning of a "vast" amount
c                of core which it can set up as it likes.
c
c
c     call symm(...)
c
c***references
c***routines called
c     a brief description of the routines in this package is given
c     below.  routines are grouped together according to the mainline
c     routine with which the're associated.  general utility routines
c     are described with symm.
c
c
c     symm  ...  the main driver routine.  initialization is done here
c                and the mainline routines ptgrp and oper are called.
c
c
c     these routines are originally the work of douglas defrees and are
c     adapted from gaussian 82.
c     ptgrp ...  determines the point group of the molecule and imposes
c                a standard orientation in cartesian space.
c
c     center  ...  determine the coordinates of the center of charge.
c     cirset  ...  search for circular sets of atoms.
c     equiv   ...  test two sets of coordinates for equivalence.
c     findc2  ...  search for a set of c2 axes perpindicular to the
c                  principal symmetry axis.
c     findcn  ...  determine the order of the principal rotation axis
c                  in a symmetric top molecule.
c     findv   ...  search for a vertical mirror plane.
c     invert  ...  inverts the molecule through the origin and return
c                  the transformation matrix.
c     or3mom  ...  calculate the third moment of charge.
c     oraxis  ...  fix the alignment of a smmetry axis with a
c                  cartesian axis.
c     orc2v   ...  orient planar c2v molecules.
c     orcn    ...  orient cs, cn, cnh, sn, and i molecules.
c     ord2h   ...  orient planar d2h molecules.
c     ordn    ...  orient dn, dnd, dnh, cnv, t, td, th, o, and oh
c                  molecules.
c     ordoc   ...  a dummy routine outlining the orientation
c                  conventions.
c     orkey   ...  determine the key atom in a symmetric top.
c     ornax   ...  determine which cartesian axis passes through the
c                  greatest number of atoms or bonds.
c     orplan  ...  orient the molecule in a cartesian plane.
c     orptst  ...  determine if a molecule is contained in a
c                  cartesian plane and which one if it is.
c     oryz    ...  put a planar molecule in the yz plane and orient
c                  it.
c     reflct  ...  reflects the molecule through one of the three cartesian
c                  planes and returns the transformation matrix.
c     rotate  ...  rotates the molecule through a given angle about one
c                  of the cartesian axes and returns the transformation
c                  matrix.
c     secmom  ...  calculate the principal second moments and axes of
c                  charge.
c     sphere  ...  determine if a spherical top molecule is
c                  tetrahedral, octahedral, or icosahedral and do som
c                  preliminary orientation.
c     sphset  ...  search for spherical sets of atoms.
c     triang  ...  given three atoms, detrmine the sides and angles
c                  of the triangle that they form.
c     tstc3   ...  do three given atoms define a c3 rotation axis?
c     tstc4   ...  do three given atoms define a c4 rotation axis?
c     tstc5   ...  do three given atoms define a c5 rotation axis?
c
c     these routines are primarily the work of paul saxe, with an
c     assist by richard martin.
c     oper  ...  generates the 3x3 transformation matrices, atom permutation
c                lists, and salc's. 
c
c     atmset  ...  finds the sets of symmetry related atoms.
c     bftoso  ...  makes the basis function to salc transformation 
c                  matrix.
c     c2prim  ...  generate the operations of a set of c2 axes int the
c                  xy plane.
c     ci      ...  generate the operation of inversion.
c     cnaxis  ...  generate the operations of a proper rotation about
c                  the z-axis (or another if overidden by "forax").
c     genrtr  ...  generates the group representation matrices and
c                  character table.
c     group   ...  determines the number of generators, the number of
c                  irreducible representations and their labels, etc.
c     permut  ...  generates the atom permutation array.
c     projct  ...  applys projection operators to a
c                  reducible representation.
c     relate  ...  determines the atoms which transform into one another
c                  under the operations of the group.
c     salc    ...  forms the symmetry-adapted-linear-combinations.
c     sigmah  ...  generate the operations produced by the product of
c                  a proper rotation and horizontal reflection.
c     sigmav  ...  generate the operations produced by the product of
c                  a proper rotation and vertical reflection.
c     snaxis  ...  generate the operations of a rotation-reflection axis.
c     sprint  ...  prints assorted information for debugging.
c     sympt   ...  generate various pointers for the integral codes.
c     symset  ...  generates pointers to the symmetry contraction sets.
c     tmat    ...  generates the transformation matrices for the 
c                  coordinates under each operation of the group.
c     tolocl  ...  rotate and translate to the local frame.
c     tomast  ...  rotate and translate to the master frame.
c     trmat   ...  generates the transformation matrices for the basis
c                  functions.
c     octa    ...  generate the operations of molecules in the point
c                  groups o and oh.
c     tetra   ...  generate the operations of molecules in the point
c                  groups t and td.
c     vert    ...  generate the operations of a set of vertical
c                  planes.
c
c***end prologue       symm
      implicit integer(a-z)
      integer ian(natoms)
      integer a
c
      character*8 symflg
      character*4 pgrp,subgrp
      character*(*) forax, forgrp
c
      logical inau, dump, norot, redsym, prtsym
c
      real*8 cin(3,natoms), atmchg(natoms), cnew(maxap3,3)
      real*8 ctmp(natoms,3)
      real*8 trvec(3), rotmat(3,3)
      real*8 one, toang, toler, tol2, told1, told2,zero
      real*8 z
c
      parameter (zero=0.0d+00,one=1.0d+00)
c
      common/io/inp,iout
      common/tol/ toler,tol2
      data told1/1.0d-08/, told2/1.0d-07/
      save told1,told2
      pointer(p,z(1)),(p,a(1))
c
 1000 format(1x,'symmetry:')
 1010 format(5x,'full point group:',2x,a4)
 1020 format(5x,'user requests   :',2x,a4)
c
c     announce to the world that we're here.
      if(prtsym) write(iout,1000)
c
c     initialize pgrp, subgrp. 
c     pgrp is a character string in which the full point group will
c     be returned.  subgrp is the largest abelian subgroup. while this
c     package will generate the salc's, transformation matrices, etc.
c     for nearly any point group, the rest of mesa can utilize only
c     information about subgroups of d2h symmetry or less.  therefore,
c     the appropriate subgroup is determined in 'group' and the 
c     information passed to the rest of the program refers to that
c     subgroup and not the full point group.
c
c     initialize toler, tol2.   
c     toler is a constant which is the
c     smallest distance in angstroms which will be considered finite.
c     tol2 is a constant which is used for comparing non-coordinate
c     quantities for zero (such as the difference between principal
c     moments of charge) and as a cutoff for the values of cartesian
c     coordinates in omega (if c(iat,ixyz) .lt. tol then c(iat,ixyz)
c     = zero).
c
      pgrp='    '
      subgrp='    '
      if (forax.eq.' ') forax='z'
      symflg='ok'
      toler = told1
      tol2  = told2
c
c     change toler if the coordinates are in atomic units.
      if(inau) toler=toler/toang
c
c     copy the input vector of cartesian coordinates to the internal
c     array.
      do 10 i = 1, natoms
         do 10 j = 1, 3
            ctmp(i,j) = cin(j,i)
   10 continue
c
c     default the new coordinates and rotation matrix.
      call rzero(cnew,3*maxap3)
      do 15 i=1,natoms
         do 15 j=1,3
            cnew(i,j)=ctmp(i,j)
   15 continue
      do 20 i=1,3
         cnew(i+natoms,i)=one
   20 continue
c
c     check for special cases.
c     is this an atomic calculation?
      if(natoms.le.1) then
         symflg='atomic'
         pgrp='o(3)'
      endif
c
c     check for dummy and ghost atoms.
      do 30 iat = 1, natoms
         if(ian(iat).le.0) symflg='ghosts'
   30 continue
c
c
c     determine the point group and standard orientation.
c     the new coordinates are returned in cnew, the schoenflies symbol
c     for the point group in pgrp, and the translation vector
c     in trvec.  the final rotation matrix is the last three "atoms"
c     in cnew.
c
c     retrieve some basis set information.
      call iosys('read integer "maximum momentum group" '//
     $           'from rwf',1,maxmom,0,' ')
      call iosys('read integer "number of basis functions" from rwf',
     $            1,nbf,0,' ')
c
c     unless instructed otherwise, determine the molecular point group. 
c  
c     don't confuse ptgrp with an atomic calculation.
      if(symflg.ne.'atomic') then
c
c        allocate core.
         pop=1
         set=pop+natoms
         scr1=iadtwp(set+natoms)
         scr2=scr1+3*maxap3
         scr3=scr2+3*maxap3
         scr4=scr3+(natoms*(natoms+1)/2)
         scr5=scr4+natoms
         scr6=scr5+natoms
         end=wpadti(scr6+natoms)
         call getmem(end,p,ngot,'symm',0)
         call ptgrp(maxap3,cnew,z(scr1),ctmp,z(scr2),ian,atmchg,
     $              natoms,dump,pgrp,trvec,symflg,z(scr3),z(scr4),
     $              z(scr5),z(scr6),a(pop),a(set))
         call getmem(-ngot,p,idum,'symm',idum)
      endif
      if(prtsym) write(iout,1010) pgrp
c
c     if the user has requested that we not use the rotation to the
c     standard orientation, restore the original coordinates and
c     default rotation matrix and translation vector.
      if(norot) then
         call rzero(cnew,3*maxap3)
         do 35 i=1,natoms
            do 35 j=1,3
               cnew(i,j)=cin(j,i)
   35    continue
         do 36 i=1,3
            trvec(i)=zero
            cnew(i+natoms,i)=one
   36    continue
      end if
c
c     if this is an atom, request d2h symmetry.
      if(symflg.eq.'atomic') then
         pgrp='d2h'
      endif
c
c     if the user has requested a specific group, override the
c     full point group.
      if (forgrp.ne.' ') then
         pgrp=forgrp
         write(iout,1020) pgrp
      end if
c
c
c     generate the permutation matrices, salc's, etc.
      call getmem(maxmom+1,p,ngot,'symm',0)
      ncart=1
      call iosys('read integer "number of cartesians" '//
     $           'from rwf',maxmom+1,a(ncart),0,' ')
      call oper(maxap3,pgrp,subgrp,forax,natoms,maxop,nop,cnew,ctmp,
     $          prtsym,dump,redsym,maxmom,a(ncart),nbf,atmchg)
      call getmem(-ngot,p,idum,'symm',idum)
      numop=nop
c
c     return to the original orientation.
      do 50 i=1,3
         do 40 j=1,3
            rotmat(i,j)=cnew(natoms+j,i)
  40     continue
  50  continue
      call deornt(rotmat,0)
c
c
      return
      end
