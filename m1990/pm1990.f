*deck @(#)pm1990.f	1.1  11/30/90
      subroutine pm1990(z,a,maxcor)
c***begin prologue     m1990
c***date written       891210     (yymmdd)
c***keywords           m1990, link 1990, graphics, integrals, property
c***author             martin, richard    (lanl)
c***source             m1990
c***purpose            computes the value of the primitive basis functions
c                      on a grid of points.
c***description
c
c***references
c
c***routines called
c                      m1990
c
c***end prologue       m1990
      implicit integer(a-z)
      real*8 z(*)
      integer a(maxcor)
      parameter (maxpts=128)
      real*8 gridx(maxpts),gridy(maxpts),gridz(maxpts)
      real*8 step,fpkey
      character*4096 ops
      character*16 chknam
c
c     a and z are equivalenced in the calling routine.
c
      common/io/inp,iout
c
c     retrieve the options string.
      call iosys('read character options from rwf',0,0,0,ops)
c
c     determine the step size for the grid; default=0.3au.
      step=fpkey(ops,'graphics=stepsize=',0.3d+00,' ')
c
c     open the checkpoint file.
      call iosys('read character "checkpoint filename" from rwf',
     $            0,0,0,chknam)
      call iosys('open chk as old',1000000,0,0,chknam)
c
c     get the lengths of arrays needed for core allocation.
c     ntypes is the total number of types used to define the lengths
c     in m102.  this total is composed of two sets: the first nbtype
c     types are used to mark the basis function types.  the remainder
c     refer to the ecp types.
      call iosys('read integer "number of atoms" from chk',1,nat,0,' ')
      call iosys('read integer "number of basis functions" from chk',1,
     $            nbasis,0,' ')
      call iosys('length of exponents on chk',nprim,0,' ')
      call iosys('read integer "number of primitive functions"'
     $            //' from chk',1,npf,0,' ')
      call iosys('length of "contraction coefficients" on chk',
     $            ncont,0,' ')
      call iosys('length of "number of pure functions" on chk',
     $           ntypes,0,' ')
      call iosys('read "number basis types" from chk',
     $           1,nbtype,0,' ')
      call iosys('length of "power of x" on chk',lenxyz,0,' ')
      nnp=(nbasis+1)*nbasis/2
c
c     divide core for basis set information.
      ian=1
      c=iadtwp(ian+nat)
      cont=c+3*nat
      ex=cont+ncont
      phix=ex+nprim
      phiy=phix+maxpts*npf
      phiz=phiy+maxpts*npf
      t1=phiz+maxpts*npf
      t2=t1+nbasis*nbasis
      mo=t2+npf*nbasis
      lowdin=mo+npf*nbasis
      s=lowdin+npf*nbasis
      smhalf=s+nnp
      u=smhalf+nnp
      eigval=u+nbasis*nbasis
      scr1=eigval+nbasis
      scr2=scr1+nbasis*nbasis
      triang=scr2+nbasis*nbasis
      ptprim=wpadti(triang+nbasis*nbasis)
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
      need=nz+lenxyz
c
c
c
      top=need+max(top1,top2)
      if (top.gt.maxcor) then
c        call getscm(top,z(1),ngot,'m1990: main',0)
c        maxcor=maxcor+ngot
         call lnkerr('not enough core in main')
      endif
c
c     read in basis set information from checkpoint file.
      call iosys('read real exponents from chk',-1,z(ex),0,' ')
      call iosys('read real "contraction coefficients" from chk',
     $            -1,z(cont),0,' ')
      call iosys('read integer "atomic numbers" from chk',
     $           -1,a(ian),0,' ')
      call iosys('read real coordinates from chk',-1,z(c),0,' ')
      call iosys('read integer "pointer to primitives" from chk',
     $             -1,a(ptprim),0,' ')
      call iosys('read integer "number of primitives" from chk',
     $             -1,a(noprim),0,' ')
      call iosys('read integer "pointer to contraction coefficients"'//
     $            'from chk',-1,a(ptcont),0,' ')
      call iosys('read integer "number of contraction coefficients"'//
     $            'from chk',-1,a(nocont),0,' ')
      call iosys('read integer "number of cartesians" from chk',
     $            -1,a(nocart),0,' ')
      call iosys('read integer "number of pure functions" from chk',
     $            -1,a(nobf),0,' ')
      call iosys('read integer "minimum momentum" from chk',-1,
     $           a(minmom),0,' ')
      call iosys('read integer "maximum momentum" from chk',-1,
     $           a(maxmom),0,' ')
      call iosys('read integer "pointer to cartesians" from chk',
     $             -1,a(mintyp),0,' ')
      call iosys('read integer "pointer to first function" from chk',
     $            -1,a(start),0,' ')
      call iosys('read integer "power of x" from chk',-1,a(nx),0,' ')
      call iosys('read integer "power of y" from chk',-1,a(ny),0,' ')
      call iosys('read integer "power of z" from chk',-1,a(nz),0,' ')
c
c     read molecular information from the chk file.
      call iosys('read integer "number of alpha electrons" from chk',
     $           1,nae,0,' ')
      call iosys('read integer "number of beta electrons" from chk',
     $           1,nbe,0,' ')
      call iosys('read integer "charge" from chk',1,charge,0,' ')
      call iosys('read integer "spin multiplicity" from chk',
     $            1,multip,0,' ')
c
c     read in the molecular orbitals.
      call iosys('read real "overlap integrals" from chk',
     $           -1,z(s),0,' ')
      call iosys('read real "scf vector" from chk',-1,z(t1),0,' ')
c
c     read in the transformation from primitives to basis functions.
      call iosys('read real t(prim,cont) from chk',-1,z(t2),0,' ')
c
c     close the checkpoint file.
      call iosys('close chk',namchk,0,0,' ')
c
c     generate the grid.
      call grid(gridx,gridy,gridz,ngridx,ngridy,ngridz,
     $          maxpts,xmin,ymin,zmin,incrx,incry,incrz,
     $          z(c),nat,step)
c
c     calculate the fermi contact integrals over the grid.
      call dirac(z(c),z(ex),z(cont),a(ptprim),a(noprim),
     #           a(nocont),a(ptcont),nat,nprim,
     #           ntypes,nbtype,nnp,ncont,a(start),nbasis,
     #           a(nocart),a(nobf),a(maxmom),a(minmom),a(mintyp),
     #           a(nx),a(ny),a(nz),ops,gridx,gridy,gridz,
     #           ngridx,ngridy,ngridz,z(phix),z(phiy),z(phiz),npf)
c
c     write this information to a formatted file for processing
c     on the titan.
c
      nmo=nbasis
c     transform the mo's to the primitive basis.
c     the transformation matrix t(prim,cont) has both the contraction
c     coefficients and the normalization constants built in(see m102).
      call ebc(z(mo),z(t2),z(t1),npf,nbasis,nmo)
c
c     find the lowdin orthogonal molecular orbitals.
c     first find s**-1/2 through the call to sinv.
      iprint=0
      call sinv(z(s),z(smhalf),z(u),z(eigval),z(scr1),z(scr2),
     $          nbasis,nnp,z(triang),iprint)
c     get s**1/2 = s*s**-1/2
      call trtosq(z(scr1),z(s),nbasis,nnp)
      call trtosq(z(scr2),z(smhalf),nbasis,nnp)
      call ebc(z(u),z(scr1),z(scr2),nbasis,nbasis,nbasis)
c     form (s**+1/2)*c in the contracted basis.
      call ebc(z(scr1),z(u),z(t1),nbasis,nbasis,nmo)
c     transform this matrix to the primitive basis.
      call ebc(z(lowdin),z(t2),z(scr1),npf,nbasis,nmo)
c
      call xyzout(z(phix),z(phiy),z(phiz),ngridx,ngridy,ngridz,nat,npf,
     $            nmo,z(mo),z(lowdin),z(c),a(ian),xmin,ymin,zmin,
     $            incrx,incry,incrz,nbasis,charge,multip,nae+nbe,nae,
     $            nbe)
c
c
      call chainx(0)
      return
      end
