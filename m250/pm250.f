*deck @(#)pm250.f	5.1  11/6/94
      subroutine pm250(z,a)
c***begin prologue     pm250.f
c***date written       880201  
c***revision date      11/6/94      
c
c***keywords           reaction path hamiltonian
c***author             page, michael(nrl)
c***source             @(#)pm250.f	5.1   11/6/94
c***purpose            
c***description
c   this program computes the steepest descent reaction
c   pathway and all of the quantities needed to construct
c   the reaction path hamiltonian of miller, handy and adams.
c
c***references
c   w. miller, n. handy and j. adams, j. chem. phys. v72 p90  (1980)
c   m. page and j. mciver,            j. chem. phys. v88 p922 (1988)
c
c***routines called
c
c***end prologue       pm250.f
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
      integer nvar,natoms,nz,nvv,nvvc
      integer atmass,cmass,ftric,kian,order,xc,gc
      integer fsq,zf,top
      integer wpadti,iadtwp
      integer maxcor
      character*4096 ops
      character*16 vname(500)
      character*4 ians
      logical rphdone
      real*8 stotal,rphstep
c
      common/io/inp,iout
c
c     --- get the link options.
      call iosys('read character options from rwf',-1,0,0,ops)
      call iosys('read integer "number of variable coordinates" '//
     $           'from rwf',1,nvar,0,' ')
      call iosys('read integer "number of z-matrix entries" from rwf',
     $            1,nz,0,' ')
      call iosys('read integer "number of atoms" from rwf',1,natoms,
     $           0,' ')
c
      call iosys('does rph_status exist on rwf',0,0,0,ians)
      if(ians.eq.'no') then
         rphdone=.false.
      else
         call iosys('read integer rph_status from rwf',1,rphdone,0,' ')
      endif
c
c     --- allocate core
      nvvc=3*natoms*(3*natoms+1)/2
      nvv=nvar*(nvar+1)/2
      atmass=1
      cmass=atmass+natoms
      ftric=cmass+3*natoms
      kian=wpadti(ftric+nvvc)
      order=kian+natoms
      xc=iadtwp(order+3*natoms)
      gc=xc+3*natoms
      fsq=gc+3*natoms
      zf=fsq+(3*natoms)*(3*natoms)
      top=wpadti(zf+10000)
      call getscm(top,z,maxcor,'m250',0)
c
c
      call rphdrv(nvar,nz,natoms,ops,vname,a(kian),z(atmass),
     $            z(cmass),z(xc),z(gc),z(zf),z(ftric),z(fsq),nvvc,
     $            a(order),rphdone,rphstep,stotal)
      call iosys('write integer rph_status to rwf',1,rphdone,0,' ')
c
c     --- if we're not finished, override the jump in the route.
c         if we are finished or have an exit flag for some other 
c         reason, execute the jump
      call iosys('write real rph_stotal to rwf',1,stotal,0,' ')
      call iosys('write integer rph_step to rwf',1,rphstep,0,' ')
      if(rphdone) then
            call chainx(0)
      else
            call chainx(1)
      endif
c
c
      stop
      end
