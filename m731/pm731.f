*deck @(#)pm731.f	5.2  4/18/95
      subroutine pm731(z,a)
c***begin prologue     pm731.f
c***date written       851001   (yymmdd)
c***revision date      4/18/95
c
c   26 january  1988   bhl at brl
c       read "cartesian first derivatives"
c
c   17 february 1987   pws at lanl
c       reworking the core allocation to handle reals and integers on 32 bit
c       machines.
c
c***keywords           m731, link 731, internal forces, gradient,
c                      wilson matrices
c***author             binkley, et al., gaussian82
c                      martin, richard    (lanl)
c***purpose
c***description
c***references
c***routines called
c***end prologue       pm731.f
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
      integer inp,iout,maxcor
      integer ian,ianz,iz,bl,alpha,beta,lbl,lalpha,lbeta
      integer fx,fi,f,b,ib,g,xm,cz,iscr,scr,top
      integer nz,natoms,nvar,nparm
      integer wpadti,iadtwp
      character*4096 ops
      character*8 prtflg
      character pgrp*4
      logical logkey,dump,prtfor,usesym
      real*8 toang,one,fmax,frms,fmaxi,frmsi,arrmax,sdot,det,cutoff
      real*8 t(3,3)
c
      parameter (one=1.0d+00,cutoff=1.0d-12)
      data prtfor/.true./
      save prtfor
c
      common/io/inp,iout
c
 1000 format(1x,'m731:transform forces from cartesian',
     $          ' to internal coordinates')
 1005 format(5x,'max ',f12.6,5x,'rms ',f12.6)
 1010 format(' m731: ***warning*** formbg det=',e20.10)
 1020 format(5x,'axes rotated back to original set.')
c
c
c     --- get the link options.
      call iosys('read character options from rwf',-1,0,0,ops)
c
c     --- set the local options.
      dump=.false.
      if(logkey(ops,'dump',.false.,' ')) dump=.true.
c
c     --- has printing been turned off externally?
      call iosys('read character "print flag" from rwf',-1,0,0,prtflg)
      if(prtflg.eq.'minimum') prtfor=.false.
      if(logkey(ops,'print=force',.false.,' ')) prtfor=.true.
c
c     --- get pertinent variables.
      call iosys('read integer "number of z-matrix entries" from rwf',
     $            1,nz,0,' ')
      call iosys('read integer "number of atoms" from rwf',
     $            1,natoms,0,' ')
      call iosys('read integer "number of variable coordinates" '//
     $           'from rwf',1,nvar,0,' ')
      call iosys('read real angstrom/bohr from rwf',1,toang,0,' ')
c
c     --- set the number of independent internal variables.
      call iosys('read character "point group" from rwf',-1,0,0,pgrp)
      if(pgrp.eq.'d*h'.or.pgrp.eq.'c*v') then
c        --- linear molecules.
         nparm=3*nz-5
      else
         nparm=3*nz-6
      endif
c
c     --- is symmetry in use?
      call iosys('read integer usesym from rwf',-1,usesym,0,' ')
c
c     --- allocate core.
      ian=1
      ianz=ian+natoms
      iz=ianz+nz
      bl=iadtwp(iz+4*nz)
      alpha=bl+nz
      beta=alpha+nz
      lbl=wpadti(beta+nz)
      lalpha=lbl+nz
      lbeta=lalpha+nz
      fx=iadtwp(lbeta+nz)
      fi=fx+3*natoms
      f=fi+nparm
      b=f+nvar
      ib=wpadti(b+3*4*nparm)
      g=iadtwp(ib+4*nparm)
      xm=g+nparm*nparm
      cz=xm+5*nz
      iscr=wpadti(cz+3*nz)
      scr=iadtwp(iscr+nparm)
      top=wpadti(scr+4*nparm)
      call getscm(top,z,maxcor,'m731',1)
c
c     --- retrieve the cartesian energy derivatives, rotate
c         them back to the original coordinate system if symmetry is
c     in use, convert them to forces, and print.
      call iosys('read real "cartesian first derivatives" from rwf',
     $           3*natoms,z(fx),0,' ')
c
      if(prtfor) write(iout,1000)
      if(usesym) then
         if(prtfor) write(iout,1020)
         call iosys('read real "rotation matrix" from rwf',-1,t,0,' ')
         call rotf(natoms,t,z(fx),z(fx))
      endif
      call iosys('read integer "atomic numbers" from rwf',
     $     natoms,a(ian),0,' ')
      call smul(z(fx),z(fx),-one,3*natoms)
      fmax=arrmax(z(fx),3*natoms)
      frms=sqrt(sdot(3*natoms,z(fx),1,z(fx),1)/float(3*natoms))
c
c     --- print the cartesian forces.
      if(prtfor) then
         call fcorpr(iout,natoms,a(ian),z(fx))
         write(iout,1005) fmax,frms
      endif
c
c     --- retrieve the z-matrix information from the rwf.
      call rzmat('rwf',nz,nvar,a(ianz),a(iz),z(bl),z(alpha),z(beta),
     $                 a(lbl),a(lalpha),a(lbeta))
c
c     --- read the full coordinate list (dummies included).
      call iosys('read real coords&dummies from rwf',3*nz,z(cz),0,' ')
c$$$ comment this out --- the call is incorrect!  Fix that someday
c$$$      if(dump) then
c$$$         call zprint(nz,a(ianz),a(iz),z(bl),z(alpha),
c$$$     $               z(beta),toang)
c$$$      endif

c
c     --- form the wilson b and g matrices.
      call formbg(nz,a(ianz),a(iz),z(bl),z(alpha),z(beta),nparm,z(b),
     $            a(ib),z(g),z(xm),z(cz),a(iscr),z(scr),dump,toang,det)
      if(det.lt.cutoff) write(iout,1010) det
c
c     --- transform cartesian forces in fx into internal forces in f.
      call tranf(nparm,nz,natoms,a(ianz),z(fx),z(fi),
     $           a(ib),z(b),z(g),a(iscr))
      fmaxi=arrmax(z(fi),nparm)
      frmsi=sqrt(sdot(nparm,z(fi),1,z(fi),1)/float(nparm))
c
c     --- print the internal forces.
      if(prtfor) then
         call fzprnt(nz,a(ianz),a(iz),z(fi))
c        write(iout,1005) fmaxi,frmsi
      endif
c
c     --- prepare forces for optimization controller.
      call putf(nz,a(lbl),a(lalpha),a(lbeta),nparm,nvar,
     $          z(fi),z(f),dump)
      call iosys('write real forces to rwf',nvar,z(f),0,' ')
c
c     --- exit gracefully.
      call chainx(0)
c
c
      stop
      end
