*deck @(#)pm1022.f	5.1  11/6/94
      subroutine pm1022(z,a)
c
c***begin prologue     m1022
c***date written       871128  (yymmdd)
c***revision date      yymmdd  (yymmdd)
c***keywords           cphf contribution to second derivatives
c***author             saxe, paul (lanl)
c***source             @(#)pm1022.f	5.1   11/6/94
c***purpose            form the cphf contribution to the hf second
c                      derivatives.
c***description
c
c***references         "unified theoretical treatment of analytic first
c                       and second derivatives in open-shell hartree-
c                       fock theory", y. osamura, y. yamaguchi, p. saxe,
c                       m. a. vincent, j. f. gaw and h. f. schaefer iii,
c                       chemical physics 72 (1983) 131-139.
c
c***routines called
c
c***end prologue       m1022
c
      implicit integer (a-z)
c
      real*8 z(*)
      integer a(*)
      character*4096 ops
c
      common /io/ inp,iout
c
c     ----- recover the options string -----
c
      call iosys('read character options from rwf',-1,0,0,ops)
c
c     ----- recover dimensions for core allocation -----
c
      call iosys('read integer "number of basis functions" from rwf',
     $           1,num,0,' ')
      call iosys('read integer "number of shells" from rwf',
     $           1,nshell,0,' ')
      call iosys('read integer "number of independent rotations" '//
     $           'from rwf',1,nindep,0,' ')
      call iosys('read integer "number of atoms" from rwf',
     $           1,natoms,0,' ')
      nder=3*natoms
      nd2e=nder*(nder+1)/2
      nnp=num*(num+1)/2
c
c     ----- allocate core -----
c
      pt=1
      minshl=pt+nnp
      maxshl=minshl+nshell
      numshl=maxshl+nshell
      d2e=iadtwp(numshl+nshell)
      ld2e=d2e+nd2e
      ds=ld2e+nd2e
      udep=ds
      u=ds+nnp*nder
      dlag=u
      sa=u+num*num*nder
      w=sa+num**2
      need=wpadti(w+num*num*nder)
c
      call getscm(need,a,maxcor,'core',0)
c
c     ----- read in number of orbitals of shell, ... -----
c
      call iosys('read integer "number in shell" from rwf',nshell,
     #            a(numshl),0,' ')
      call iosys('read integer "first of shell" from rwf',nshell,
     #            a(minshl),0,' ')
      call iosys('read integer "last of shell" from rwf',nshell,
     #            a(maxshl),0,' ')
      call iosys('read integer "independent rotations" from rwf',nnp,
     #            a(pt),0,' ')
c
c     ----- and fill out the second derivatives -----
c
      call hfd2e(z(ds),nnp,nder,z(udep),nindep,z(u),num,z(w),z(dlag),
     $     z(sa),z(d2e),nd2e,z(ld2e),a(pt),a(minshl),a(maxshl),
     $     a(numshl),nshell,ops)
c
c     ----- and exit this link -----
c
      call chainx(0)
c
c
      stop
      end
