*deck @(#)pm732.f	5.1  11/6/94
      subroutine pm732(z,a)
c***begin prologue      pm732.f
c***date written        880114   (yymmdd)
c***revision date       11/6/94   
c***keywords
c                       m732, link 732 transform force constants
c***author              page, michael (nrl)
c***source              @(#)pm732.f	5.1 11/6/94
c***purpose
c***                    transform cartesian force constant matrix to
c***                    internal z-matrix coordinates
c***description
c
c***references
c
c***routines called	
c***end prologue 	pm732.f
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
      integer inp,iout
      integer nvar,nz,natoms
      integer nvv,nvvc,x,f,frcnst,lbl,xx,c,fsqc,ftric,fsqi
      integer kian,katchg,kc,kianz,kiz,kbl,kalpha,kbeta,klbl
      integer klalph,klbeta
      integer iscr1,iscr2,iscr3,iscr4,iscr5,iscr6,iscr7,iscr8,iscr9
      integer iscr10,iscr11
      integer bplus,bminus,bder,tt,gc,top
      integer wpadti,iadtwp,maxcor
      character*4096 ops
      character*16 vname(500)
      real*8 toang
c
      common/io/inp,iout
c
 1000 format(' in m732 : top,maxcor ',2i20)
c
c     --- get the link options.
      call iosys('read character options from rwf',-1,0,0,ops)
      call iosys('read integer "number of variable coordinates" '//
     $           'from rwf',1,nvar,0,' ')
      call iosys('read integer "number of z-matrix entries" from rwf',
     $           1,nz,0,' ')
      call iosys('read real angstrom/bohr from rwf',1,toang,0,' ')
      call iosys('read integer "number of atoms" from rwf',
     $            1,natoms,0,' ')
c
c     --- allocate core.
      nvv=nvar*(nvar+1)/2
      nvvc=3*natoms*(3*natoms+1)/2
      x=1
      f=x+nvar
      frcnst=f+nvar
      lbl=wpadti(frcnst+nvv)
      xx=iadtwp(lbl+nz)
      c=xx+nvar*(2*nvar+1)
      fsqc=c+3*natoms
      ftric=fsqc+(3*natoms)*(3*natoms)
      fsqi=ftric+nvvc
c                    integer atomic numbers(natoms).
      kian=wpadti(fsqi+nvar*nvar)
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
c                    scratch(8*nz)
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
      iscr11=iscr10+3*natoms
      bplus=iscr11+3*natoms*nvar
      bminus=bplus+nvar*(3*natoms)
      bder=bminus+nvar*(3*natoms)
      tt=bder+nvar*(3*natoms)
      gc=tt*nvar
      top=wpadti(gc+3*natoms)
c
      call getscm(top,z,maxcor,'m732',0)
      if(top.ge.maxcor) then
         write(iout,1000) top,maxcor
      endif
c
      call trmain(nvar,nvv,nz,natoms,nvvc,toang,ops,z(x),z(f),
     $            z(c),z(fsqc),z(ftric),z(xx),z(frcnst),vname,
     $            a(lbl),a(kian),a(kianz),a(kiz),z(kbl),z(kalpha),
     $            z(kbeta),a(klbl),a(klalph),a(klbeta),z(iscr1),
     $            z(iscr2),z(iscr3),z(iscr4),z(iscr5),z(iscr6),
     $            z(bplus),z(bminus),z(bder),z(tt),z(gc),
     $            z(iscr7),z(iscr8),z(iscr9),
     $            z(iscr10),z(iscr11),z(fsqi),z(katchg))
c
c
      call chainx(0)
c
c
      return
      end
