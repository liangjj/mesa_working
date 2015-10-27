*deck @(#)m205.f	5.1  11/6/94
      program m205
c***begin prologue      m205.f
c***date written        870215   (yymmdd)
c***revision date       11/6/94
c***keywords
c   m205, link 205, finite difference derivatives
c***author              page, michael (nrl)
c***purpose
c***description
c   m205 currently recognizes the option strings:
c
c   d2efd            ** calculate second derivatives of the energy
c                    ** with respect to nuclear coordinates by
c                    ** finite difference of analytical gradients
c
c   double_point     ** use two point difference formula
c                    ** default is single point
c
c   freq             ** calculate vibrational frequencies
c                    **   if the following restrictions cannot be met,
c                    **   'cartesian' should be specified:
c                    **     1. 3n-6 internal coordinates must be
c                    **        explicitly input in the z-matrix
c                    **     2. none of the first three centers can
c                    **        be dummy atoms
c
c  cartesian         ** use cartesian displacements for finite difference
c                    ** note that the geometry is still input in z-matrix
c                    ** format
c                    ** this option should be used to calculate frequencies
c                    ** if the z-matrix restrictions on 'freq' cannot be met
c
c   stpsize=a        ** stepsize for difference formula in
c                    ** radians and au. default is 0.01
c
c***references
c
c***routines called
c***end prologue 	m205.f
      implicit none
c     --- input variables -----
c     --- input arrays (unmodified) ---
c     --- input arrays (scratch) ---
c     --- output arrays ---
c     --- output variables ---
c     --- scratch arrays ---
      integer a(*)
      real*8 z(*)
      pointer (p,z(1)), (p,a(1))
c     --- local variables ---
      integer inp,iout
      integer nvv,nvvc,maxpt,x,f,xx,ff,fcc,frcnst
      integer fs,atmass,cmass,scr,scrsq,lbl,kian,katchg,kc
      integer kianz,kiz,kbl,kalpha,kbeta,klbl,klalph,klbeta
      integer iscr1,iscr2,iscr3,iscr4,iscr5,iscr6,iscr7
      integer iscr8,iscr9,iscr10,iscr11
      integer xc,fc,xxc,ffc,frcnsc,sqc,zf,top
      integer nvar,nz,natoms,wpadti,iadtwp,maxcor
      character*4096 ops
      character*16 vname(500)
      character*4 ians
      logical d2edone,abnrml
      real*8 toang
c
      data d2edone/.false./, abnrml/.false./
      save d2edone,abnrml
c
      common/io/inp,iout
      call drum
c
c     --- i don't know what this does.
      call iosys('reset alignment of rwf',0,0,0,' ')
c
c     --- get the link options.
      call iosys('read character options from rwf',-1,0,0,ops)
      call iosys('read integer "number of variable coordinates" '//
     $           'from rwf',1,nvar,0,' ')
      call iosys('read integer "number of z-matrix entries" from rwf',
     $            1,nz,0,' ')
      call iosys('read real angstrom/bohr from rwf',1,toang,0,' ')
      call iosys('read real angstrom/bohr from rwf',1,toang,0,' ')
      call iosys('read integer "number of atoms" from rwf',1,natoms,
     $            0,' ')
      call iosys('does d2e_status exist on rwf',0,0,0,ians)
      if (ians.eq.'no') then
         d2edone=.false.
      else
         call iosys('read integer d2e_status from rwf',1,d2edone,0,' ')
      end if
c
c     --- allocate core.
      nvv=nvar*(nvar+1)/2
      maxpt=6*natoms+1
      x=1
      f=x+nvar
      xx=f+nvar
      ff=xx+nvar*maxpt
      fcc=ff+nvar*maxpt
      frcnst=fcc+2*nvv
      fs=frcnst+nvv
      atmass=fs+maxpt
      cmass=atmass+natoms
      scr=cmass+3*natoms
      scrsq=scr+max(nvar,3*natoms)
      lbl=wpadti(scrsq+2*nvar*nvar)
c                    integer atomic numbers(natoms).
      kian=lbl+nz
c                    atomic charges(natoms).
      katchg=iadtwp(kian+max(natoms,nz))
c                    coordinates(natoms).
      kc=katchg+max(natoms,nz)
c                    z-matrix atomic numbers(nz).
      kianz=wpadti(kc+3*max(natoms,nz))
c                    integer components of the z-matrix(4*nz).
      kiz=kianz+nz
c                    z-matrix bond lengths(nz).
      kbl=iadtwp(kiz+4*nz)
c                    z-matrix bond angles(nz).
      kalpha=kbl+nz
c                    z-matrix dihedral angles(nz).
      kbeta=kalpha+nz
c                    z-matrix bond length map(nz).
      klbl=wpadti(kbeta+nz)
c                    z-matrix angle map(nz).
      klalph=klbl+nz
c                    z-matrix dihedral angle map(nz).
      klbeta=klalph+nz
c                    scratch
      iscr1=iadtwp(klbeta+nz)
      iscr2=iscr1+nz
      iscr3=iscr2+nz
      iscr4=iscr3+nz
      iscr5=iscr4+nz
      iscr6=iscr5+nz
      iscr7=iscr6+3*nz
      iscr8=iscr7+3*natoms*nvar
      iscr9=iscr8+3*natoms
      iscr10=iscr9+3*natoms
      iscr11=wpadti(iscr10+3*natoms)
      xc=iadtwp(iscr11+3*natoms)
      fc=xc+3*natoms
      xxc=fc+3*natoms
      ffc=xxc+3*natoms*maxpt
      frcnsc=ffc+3*natoms*maxpt
      sqc=frcnsc+3*natoms*(3*natoms+1)/2
      zf=sqc+3*natoms*3*natoms
      top=wpadti(zf+10000)
      call getmem(top,p,ngot,'m205',0)
      call izero(a,top)
c
      nvvc=3*natoms*(3*natoms+1)/2
      call d2emain(nvar,nvv,nz,natoms,maxpt,toang,ops,z(x),z(f),
     $             z(xx),z(ff),z(frcnst),z(fs),vname,a(lbl),
     $             a(kian),z(atmass),d2edone,abnrml,z(scr),
     $             z(scrsq),z(cmass),z(katchg),
     $             z(kc),a(kianz),a(kiz),z(kbl),z(kalpha),
     $             z(kbeta),a(klbl),a(klalph),a(klbeta),z(iscr1),
     $             z(iscr2),z(iscr3),z(iscr4),z(iscr5),z(iscr6),
     $             z(iscr7),z(iscr8),z(iscr9),z(iscr10),
     $             a(iscr11),z(xc),z(fc),z(xxc),z(ffc),z(zf),
     $             z(frcnsc),z(sqc),nvvc)
c
c
      call iosys('write integer d2e_status to rwf',1,d2edone,0,' ')
c
c     --- if we're not finished override the jump in the route.
c         if we have finished or have an exit flag for some other reason
c         execute the jump.
      if(d2edone.or.abnrml) then
         call chainx(0)
      else
         call chainx(1)
      endif
c
c
      stop
      end
