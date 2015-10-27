*deck @(#)pm1021.f	5.1  11/6/94
      subroutine pm1021(z,a)
c
c***begin prologue     m1021
c***date written       861211  (yymmdd)
c***revision date      yymmdd  (yymmdd)
c***keywords           u vectors, coupled perturbed hartree-fock, cphf
c***author             saxe, paul (lanl)
c***source             @(#)pm1021.f	5.1   11/6/94
c***purpose            fill out the dependent blocks of u.
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
c***end prologue       m1021
c
      implicit integer (a-z)
c
      character*4096 ops
      character*128 tints,dints,rdints
      logical mcscf,logkey
      integer a(*)
      real*8 z(*)
c
c     ----- lenval is the maximum size of the buffer to read mo integrals
c
      parameter (lenval=150000)
c
      common /io/ inp,iout
c
c     ----- recover the options string -----
c
c
      call iosys('read character options from rwf',-1,0,0,ops)
c
c     ----- recover dimensions for core allocation -----
c
      mcscf=logkey(ops,'mcscf',.false.,' ')
c
      call iosys('read integer "number of basis functions" from rwf',
     $            1,num,0,' ')
      call iosys('read integer "number of atoms" from rwf',
     $            1,natoms,0,' ')
      nder=3*natoms
      nnp=num*(num+1)/2
c
c
c     ----- open the integral units -----
c
      call iosys('read character "transformed integral filename"'
     $         //' from rwf',0,0,0,tints)
      call iosys('open tints as old',0,0,0,tints)
c
      call iosys('read character "derivative integral filename"'
     $         //' from rwf', 0,0,0,dints)
      call iosys('open dints as old',0,0,0,dints)
c
      call iosys('read character "raw derivative integral filename"'
     $         //' from rwf', 0,0,0,rdints)
      call iosys('open rdints as old',0,0,0,rdints)
c
      if(.not.mcscf) then
         call iosys('read integer "number of shells" from rwf',
     $               1,nshell,0,' ')
         call iosys('read integer "number of independent rotations" '//
     $              'from rwf',1,nindep,0,' ')
c
c        ----- allocate core for cphf -----
c
         orbshl=1
         pt=orbshl+num
         minshl=pt+nnp
         maxshl=minshl+nshell
         numshl=maxshl+nshell
         f=iadtwp(numshl+nshell)
         isq=f+nshell
c
         top=wpadti(isq)
         call getscm(top,a,maxcor,'core',0)
c
c        ----- read in number of orbitals of shell, f's,... -----
c
         call iosys('read integer "number in shell" from rwf',nshell,
     #               a(numshl),0,' ')
         call iosys('read integer "first of shell" from rwf',nshell,
     #               a(minshl),0,' ')
         call iosys('read integer "last of shell" from rwf',nshell,
     #               a(maxshl),0,' ')
         call iosys('read real f from rwf',nshell,z(f),0,' ')
         call iosys('read integer "independent rotations" from rwf',nnp,
     #               a(pt),0,' ')
c
         nco=a(numshl)
         nconum=nco*num
c
         ds=isq+num**2
         u=ds+nnp*nder
         values=u+num*num*nder
c
c        ----- work out how many triangles of mo integrals we will hold
c
         ntriang=min(nnp**2+nconum*nnp,lenval)/(nnp+nconum)
c.debug
c        ntriang=2
c        ntriang=nnp
c.debug
         if (ntriang.lt.1) 
     $      call lnkerr('cant hold one triangle of mo ints')
c
c        ----- overlap space for the lagrangians with 'values' -----
c              as well as vectors and temporary arrays
c
         intj=values+nnp*ntriang
         clag=intj+nconum*ntriang
c
         c=clag+num
         t1=c+num**2
         t2=t1+num**2
         dclag=t2+num**2
         fj=dclag+num**2*nder
         fk=fj+nnp*nder
         uk=fk+num*num*nder
         top=wpadti(uk+num*nco*nder)
c.debug
c        ints=uk+num*nco*nder
c        top=ints+num*num*num*num
c.debug
c
c
         call getscm(top,a,maxcor,'core',0)
c
c        ----- and fill out the u vectors -----
c
c        call dumbu(z(values),nnp,num,ntriang,z(isq),z(ds),z(u),z(c),
c     #             z(dclag),a(pt),nindep,z(f),nshell,a(orbshl),
c     #             a(minshl),a(maxshl),a(numshl),z(clag),
c     #             z(t1),z(t2),nder,z(ds),ops,
c     #             nco,z(fj),z(fk),z(intj),z(uk),z(ints))
c
         call fillu(z(values),nnp,num,ntriang,z(isq),z(ds),z(u),z(c),
     #              z(dclag),a(pt),nindep,z(f),nshell,a(orbshl),
     #              a(minshl),a(maxshl),a(numshl),z(clag),
     #              z(t1),z(t2),nder,z(ds),ops,
     #              nco,z(fj),z(fk),z(intj),z(uk),nconum)
c
      else
c
c        ----- allocate core for cpmcscf -----
c
         call iosys('read integer mc_ncore from rwf',1,nco,0,' ')
         call iosys('read integer mc_nactive from rwf',1,nao,0,' ')
c
         nconum=nco*num
c
         isq=1
         ds=isq+num**2
         u=ds+nnp*nder
         values=u+num*num*nder
c
c        ----- work out how many triangles of mo integrals we will hold
c
         ntriang=min(nnp**2+nconum*nnp,lenval)/(nnp+nconum)
c.debug
c        ntriang=2
c        ntriang=nnp
c.debug
         if (ntriang.lt.1) 
     $      call lnkerr('cant hold one triangle of mo ints')
c
c        ----- overlap space for the lagrangians with 'values' -----
c              as well as vectors and temporary arrays
c
         intj=values+nnp*ntriang
         clag=intj+nconum*ntriang
c
         c=clag+num
         t1=c+num**2
         t2=t1+num**2
         dclag=t2+num**2
         fj=dclag+num**2*nder
         fk=fj+nnp*nder
         uk=fk+num*num*nder
         top=wpadti(uk+num*nco*nder)
c.debug
c        ints=uk+num*nco*nder
c        top=ints+num*num*num*num
c.debug
c
c
         call getscm(top,a,maxcor,'core',0)
c
         call mcscfu(z(values),nnp,num,ntriang,z(isq),z(ds),z(u),z(c),
     #               z(dclag),z(clag),z(t1),z(t2),nder,ops,
     #               nco,z(fj),z(fk),z(intj),z(uk),nconum,nao)
c
      end if
c
c
c     ----- and exit this link -----
c
      call iosys('close tints',0,0,0,' ')
      call iosys('close rdints',0,0,0,' ')
      call iosys('close dints',0,0,0,' ')
      call chainx(0)
c
c
      stop
      end
