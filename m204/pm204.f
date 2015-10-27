*deck @(#)pm204.f	5.1  11/6/94
      subroutine pm204(z,a)
c***begin prologue	pm204.f 
c***date written        880125   (yymmdd)
c***revision date       11/6/94
c***keywords
c   m204, link 204, vibrational frequencies
c***author  		page, michael (nrl)
c***purpose
c***description
c
c
c***references
c
c***routines called
c***end prologue	pm204.f 
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
      integer nvv,atmass,cmass,ftric,kian,order,xc,fsq,zf,top
      integer nvvc
      integer maxcor
      integer wpadti,iadtwp
      character*4096 ops
      character*16 vname(500)
c
      common/io/inp,iout
c
 1000 format(' in m204 : top,maxcor',2i20)
c
c     --- get the link options.
      call iosys('read character options from rwf',-1,0,0,ops)
      call iosys('read integer "number of variable coordinates" '//
     $           'from rwf',1,nvar,0,' ')
      call iosys('read integer "number of z-matrix entries" from rwf',
     $           1,nz,0,' ')
      call iosys('read integer "number of atoms" from rwf',
     $           1,natoms,0,' ')
c
c     --- allocate core.
      nvvc=3*natoms*(3*natoms+1)/2
      nvv=nvar*(nvar+1)/2
      atmass=1
      cmass=atmass+natoms
      ftric=cmass+3*natoms
c     integer atomic numbers(natoms).
      kian=wpadti(ftric+nvvc)
      order=kian+natoms
      xc=iadtwp(order+3*natoms)
      fsq=xc+3*natoms
      zf=fsq+(3*natoms)*(3*natoms)
      top=wpadti(zf+20000)
c
      call getscm(top,z,maxcor,'m204',0)
      if(top.ge.maxcor) then
         write(iout,1000) top,maxcor
         call lnkerr('need more core in m204')
      endif
c
c
      call frmain(nvar,nz,natoms,ops,vname,a(kian),z(atmass),
     $            z(cmass),z(xc),z(zf),z(ftric),z(fsq),nvvc,
     $            a(order))
c
c
      call chainx(0)
c
c
      return
      end
