*deck @(#)pm1020.f	5.1  11/6/94
      subroutine pm1020(z,a)
c
c***begin prologue     m1020
c***date written       861213(yymmdd)
c***revision date      930711  (yymmdd)
c   july 11, 1993      rlm at lanl
c      opening the tints file as opposed to the ints file.
c***keywords           cphf equations
c***author             saxe, paul (lanl)
c***source             @(#)pm1020.f	5.1   11/6/94
c***purpose            solving the cphf equations
c***description
c
c***references
c***routines called
c***end prologue       m1020
c
      implicit integer (a-z)
c
      real*8 z(*)
      character*4096 ops
      character*128 tints
      integer a(*)
      logical logkey
c
      common /io/ inp,iout
c
c     ----- recover the options string -----
      call iosys('read character options from rwf',-1,0,0,ops)
c
c     ----- recover dimensions for core allocation -----
      call iosys('read integer "number of independent rotations" '//
     $     'from rwf',1,nindep,0,' ')
      call iosys('read integer "number of atoms" from rwf',
     $     1,natoms,0,' ')
      nder=3*natoms
c
c     ----- open the tints unit, which contains the hessian and right-
c           hand sides
      call iosys('read character "transformed integral filename"'
     $           //' from rwf',0,0,0,tints)
      call iosys('open tints as old',0,0,0,tints)
c
c     ----- decide on direct or iterative solution of the equations --
c
c
c     ----- allocate core -----
c
      iwork=1
      hessn=iadtwp(iwork+nindep)
      b0=hessn+nindep**2
      work=b0+nindep*nder
      top=wpadti(work+nindep)
c
      call getscm(0,a,maxcor,'?',0)
c
      if (top.le.maxcor.and.
     $    logkey(ops,'cphf=direct-solution',.true.,' ')) then
c
         call getscm(top,a,maxcor,'direct',0)
c
         call direct(a(iwork),z(hessn),z(b0),z(work),nindep,nder,ops)
c
      else
c
c        ----- iterative solution of the cphf equations -----
c
      end if
c
c
      call chainx(0)
c
c
      stop
      end
