*deck @(#)pm619.f	5.1 11/6/94
      subroutine pm619(z,a)
c***begin prologue     pm619.f
c***date written       860826   (yymmdd)
c***revision date      11/6/94
c   february 12, 1994  rlm at lanl
c      extracted from the properties link(m1902) a piece to generate 
c      just the electrostatic potentials on a grid. 
c***keywords           m619, link 619, electrostatic, potential, integrals,
c***author             saxe, paul and martin, richard    (lanl)
c***source             @(#)pm619.f	5.1 11/6/94
c***purpose            computes electrostatic potential integrals
c***description
c
c***references
c
c***routines called
c***end prologue       pm619.f
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
      integer maxnbf
      parameter (maxnbf=2000)
      character*4096 ops
      character*80 card
      character*16 bflabl(maxnbf)
      logical logkey,positn,grid,inau,dolp
      integer nat,nbasis,nprim,ncont,ntypes,nbtype,lenxyz,cgrid,ngrid
      integer c,ptprim,maxprm,mxcont,maxl,maxblk,npint
      integer v,vnuc,top,prmint,need,need1,need2,ngot,left
      integer ex,cont,zan,noprim,ptcont,nocont,nocart,nobf,minmom,maxmom
      integer mintyp,start,nx,ny,nz,nnp
      integer wpadti,iadtwp
      integer inp,iout
      real*8 toang
c
      common/io/inp,iout
c
c     --- collect the options string.
      call iosys('read character options from rwf',-1,0,0,ops)
c
c     --- get the basis function labels.
      call iosys('read character "basis function labels" from rwf',
     $          -1,0,0,bflabl)
c
c     --- get the lengths of arrays needed for core allocation.
c        ntypes is the total number of types used to define the lengths
c        in m102.  this total is composed of two sets: the first nbtype
c        types are used to mark the basis function types.  the remainder
c        refer to the ecp types.
      call iosys('read integer "number of atoms" from rwf',1,nat,0,' ')
      call iosys('read integer "number of basis functions" from rwf',
     $     1,nbasis,0,' ')
      call iosys('length of exponents on rwf',nprim,0,0,' ')
      call iosys('length of "contraction coefficients" on rwf',
     $     ncont,0,0,' ')
      call iosys('length of "number of pure functions" on rwf',
     $     ntypes,0,0,' ')
      call iosys('read integer "number basis types" from rwf',
     $     1,nbtype,0,' ')
      call iosys('length of "power of x" on rwf',lenxyz,0,0,' ')
      nnp=(nbasis+1)*nbasis/2
c
c     determine the room need to hold the electrostatic potential,etc.
c     if there is a $grid section in the input stream, then read the 
c     coordinates of the points, else evaluate the potential, etc.
c     at the nuclear centers. 
c
c     --- allocate some core
c         test points in z(cgrid)
      cgrid=1
      grid=positn('$grid',card,inp)
      if(grid) then
         inau=logkey(ops,'geom=inau',.false.,' ')
         call iosys('read real angstrom/bohr from rwf',1,toang,0,' ')
         call gcoord(ngrid,z(cgrid),inau,toang) 
      else if(logkey(ops,'cox-williams',.false.,' ')) then
         call iosys('read integer "number of solvent surface points"'
     $            //' from rwf',1,ngrid,0,' ')
         call iosys('read real "solvent surface coordinates" from rwf',
     $              3*ngrid,z(cgrid),0,' ')
      else
c        evaluate at the nuclear centers. note that ngrid is
c        set here so it does not include any point charges in
c        the atoms list.
         call iosys('read integer '//
     $              '"number of atoms with basis functions"'//
     $              ' from rwf',1,ngrid,0,' ')
c        --- copy atomic coordinates into cgrid array.
         call scopy(3*ngrid,z(c),1,z(cgrid),1)
      endif
c
c     --- divide core for basis set information.
      zan=cgrid+3*ngrid
      c=zan+nat
      cont=c+3*nat
      ex=cont+ncont
      ptprim=wpadti(ex+nprim)
      noprim=ptprim+ntypes*nat
      nocont=noprim+ntypes*nat
      ptcont=nocont+ntypes*nat
      start=ptcont+ntypes*nat
      nocart=start+ntypes*nat
      nobf=nocart+ntypes
      minmom=nobf+ntypes
      maxmom=minmom+ntypes
      mintyp=maxmom+ntypes
      nx=mintyp+ntypes
      ny=nx+lenxyz
      nz=ny+lenxyz
      top=nz+lenxyz
c
c        --- retrieve information about the most demanding shell block.
      call iosys('read integer maxprm from rwf',1,maxprm,0,' ')
      call iosys('read integer maxcont from rwf',1,mxcont,0,' ')
      call iosys('read integer maxl from rwf',1,maxl,0,' ')
      call iosys('read integer maxblk from rwf',1,maxblk,0,' ')
      call iosys('read integer dolp from rwf',1,dolp,0,' ')
c
c     --- electrostatic potential integrals need scratch
      vnuc=iadtwp(top)
      v=vnuc+ngrid
      top=wpadti(v+nnp*ngrid)
c
c     --- determine the amount of scratch space required.
      npint=maxprm*maxprm
      prmint=1+npint*(16)+npint*(maxl+maxl+1)
      need1=wpadti(prmint+npint*(7*maxblk+9*(maxl+1)*(maxl+1)+1))
      need1=max(need1,wpadti(nbasis*nbasis))
      need2=wpadti(prmint+maxblk*(7*npint+mxcont*mxcont) +mxcont*maxprm)
c
      need=top+max(need1,need2)
      call getscm(need,z(1),ngot,'m619: main',0)
      left=iadtwp(ngot-need)
c
c     --- read in basis set information from read-write file.
      call iosys('read real exponents from rwf',-1,z(ex),0,' ')
      call iosys('read real "contraction coefficients" from rwf',
     $     -1,z(cont),0,' ')
      call iosys('read real "nuclear charges" from rwf',
     $           -1,z(zan),0,' ')
      call iosys('read real coordinates from rwf',-1,z(c),0,' ')
      call iosys('read integer "pointer to primitives" from rwf',
     $     -1,a(ptprim),0,' ')
      call iosys('read integer "number of primitives" from rwf',
     $     -1,a(noprim),0,' ')
      call iosys('read integer "pointer to contraction coefficients"'//
     $     ' from rwf',-1,a(ptcont),0,' ')
      call iosys('read integer "number of contraction coefficients" '//
     $     'from rwf',-1,a(nocont),0,' ')
      call iosys('read integer "number of cartesians" from rwf',
     $     -1,a(nocart),0,' ')
      call iosys('read integer "number of pure functions" from rwf',
     $     -1,a(nobf),0,' ')
      call iosys('read integer "minimum momentum" from rwf',
     $     -1,a(minmom),0,' ')
      call iosys('read integer "maximum momentum" from rwf',
     $     -1,a(maxmom),0,' ')
      call iosys('read integer "pointer to cartesians" from rwf',
     $     -1,a(mintyp),0,' ')
      call iosys('read integer "pointer to first function" from rwf',
     $     -1,a(start),0,' ')
      call iosys('read integer "power of x" from rwf',-1,a(nx),0,' ')
      call iosys('read integer "power of y" from rwf',-1,a(ny),0,' ')
      call iosys('read integer "power of z" from rwf',-1,a(nz),0,' ')
c
c
      call v0(z(c),z(ex),z(top),z(cont),a(ptprim),a(noprim),
     $       a(nocont),a(ptcont),nat,nprim,left,
     $       ntypes,nbtype,nnp,ncont,a(start),nbasis,z(zan),
     $       a(nocart),a(nobf),a(maxmom),a(mintyp),a(nx),a(ny),
     $       a(nz),a(minmom),ops,z(v),z(vnuc),dolp,bflabl,
     $       ngrid,z(cgrid))
c
c     --- write the potentials to the rwf ---
      call iosys('write real "solvent surface electrostatic potential'
     $         //' matrices" to rwf',ngrid*nnp,z(v),0,' ')
      call iosys('write real "solvent surface nuclear electrostatic'
     $        //' potentials" to rwf',ngrid,z(vnuc),0,' ')
c
c     --- exit gracefully.
      call chainx(0)
c
c
      return
      end
