*deck @(#)pm1012.f	5.1  11/6/94
      subroutine pm1012(z,a)
c***begin prologue     m1012
c***date written       861210  (yymmdd)
c***revision date      930711  (yymmdd)
c   july 11, 1993      rlm at lanl
c      removing the open of the ints file.
c***keywords           hessian, coupled perturbed hartree-fock, cphf
c***author             saxe, paul (lanl)
c***source             @(#)pm1012.f	5.1   11/6/94
c***purpose            form the hf hessian matrix for cphf equations.
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
c***end prologue       m1012
c
      implicit integer (a-z)
c
      real*8 z(*)
      character*4096 ops
      character*128 tints
      logical logkey
      integer a(*)
c
c     ----- lenbin is the length of the arrays passing elements to sorter
c
      parameter (lenbin=1024)
c
c     ----- lenval is the maximum size of the buffer to read mo integrals
c
      parameter (lenval=50000)
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
     $     1,num,0,' ')
      call iosys('read integer "number of shells" from rwf',
     $     1,nshell,0,' ')
      nnp=num*(num+1)/2
c
c     ----- allocate core -----
c
      orbshl=1
      lab=orbshl+num
      bin=lab+lenbin
      pt=bin+lenbin
      minshl=pt+nnp
      maxshl=minshl+nshell
      numshl=maxshl+nshell
      alpha=iadtwp(numshl+nshell)
      beta=alpha+nshell**2
      val=beta+nshell**2
      isq=val+lenbin
      values=isq+num**2
c
c     ----- work out how many triangles of mo integrals we will hold
c
      ntriang=min(nnp**2,lenval)/nnp
      if (ntriang.lt.1) call lnkerr('cant hold one tirangle of mo ints')
c
      core=values+nnp*ntriang
c
c     ----- overlap space for the lagrangians with values -----
c
      lag=values
      glag=lag+num**2
      top=wpadti(core+nnp**2)
c
c     ----- find the amount of core available -----
      call getscm(0,z,maxcor,'maximum available')
      if(top.gt.maxcor) then
         call lnkerr('m1012 need more core')
      endif
c
      icore=wpadti(core)
      lencor=maxcor-icore
c
c     ----- read in number of orbitals of shells, alphas and betas ----
c
      call iosys('read integer "number in shell" from rwf',nshell,
     #            a(numshl),0,' ')
      call iosys('read integer "first of shell" from rwf',nshell,
     #            a(minshl),0,' ')
      call iosys('read integer "last of shell" from rwf',nshell,
     #            a(maxshl),0,' ')
      call iosys('read real alpha from rwf',nshell**2,z(alpha),0,' ')
      call iosys('read real beta from rwf',nshell**2,z(beta),0,' ')
c
c     ----- create the array of pointers to the independent elements
c
      call ptrs(a(minshl),a(maxshl),a(numshl),nshell,a(pt),nnp,nindep,
     #          a(orbshl),num)
c
c     ----- open the integral units -----
c
      call iosys('read character "transformed integral filename"'
     $         //' from rwf',0,0,0,tints)
      call iosys('open tints as old',0,0,0,tints)
c
c     ----- and make the hessian -----
c
      call hessn(z(values),nnp,num,ntriang,z(isq),a(lab),z(val),a(bin),
     #           lenbin,a(pt),nindep,z(alpha),z(beta),nshell,z(core),
     #           a(icore),lencor,a(orbshl),a(minshl),a(maxshl),
     #           a(numshl),z(lag),z(glag))
c
c     ----- temporarily print the hessian -----
c
      if (logkey(ops,'print=cphf=hessian',.false.,' ')) then
         need=wpadti(1+nindep**2)
         if (need.le.maxcor) then
            call iosys('read real "scf hessian" from tints',
     $           nindep**2,z,0,' ')
            write (iout,20)
 20         format (/,10x,'the cphf hessian matrix:')
            call matout(z,nindep,nindep,nindep,nindep,iout)
         end if
      end if
c
c
c     ----- and exit this link -----
c
      call iosys('close tints',0,0,0,' ')
      call chainx(0)
c
c
      stop
      end
