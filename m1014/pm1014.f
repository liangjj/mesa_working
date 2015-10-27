*deck @(#)pm1014.f	5.1  11/6/94
      subroutine pm1014(z,a)
c***begin prologue     m1014
c***date written       861211  (yymmdd)
c***revision date      yymmdd  (yymmdd)
c***keywords           second derivatives, omega terms
c***author             saxe, paul (lanl)
c***source             @(#)pm1014.f	5.1   11/6/94
c***purpose            form the omega(a's) [wa's] for hf second
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
c***end prologue       m1014
c
      implicit integer (a-z)
c
      character*4096 ops
      character*128 tints
      logical logkey
      integer a(*)
      real*8 z(*)
c
c     ----- lenval is the maximum size of the buffer to read mo integrals
c
c..temp      parameter (lenval=30000)
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
      call iosys('read integer "number of independent rotations" '//
     $     'from rwf',1,nindep,0,' ')
      call iosys('read integer "number of atoms" from rwf',
     $     1,natoms,0,' ')
      nder=3*natoms
      nnp=num*(num+1)/2
c
c     ----- allocate core -----
c
      orbshl=1
      pt=orbshl+num
      minshl=pt+nnp
      maxshl=minshl+nshell
      numshl=maxshl+nshell
      alpha=iadtwp(numshl+nshell)
      beta=alpha+nshell**2
      isq=beta+nshell**2
      ds=isq+num**2
c     wa=ds+num*num*nder
      wa=ds+nnp*nder
      values=wa+num*num*nder
c
c     ----- work out how many triangles of mo integrals we will hold
c
c..temp.dumb
      if(logkey(ops,'m1014=debug',.false.,' ')) then
       lenval=nnp*nnp
      else
       lenval=max(nnp,40000)
      end if
c..temp.dumb
      ntriang=min(nnp**2,lenval)/nnp
c..temp
      if(ntriang.lt.1)
     # call lnkerr('cant hold 1 triangle of mo ints')
c..temp
c
c     ----- overlap space for the lagrangians with 'values' -----
c            as well as vectors and temporary arrays
c
      lag=values
      glag=lag+num**2
c
      t1=values
      t2=t1+num**2
      dlag=t2+num**2
c
      top=wpadti(max(values+nnp*ntriang,glag+num**2*nshell,
     +               dlag+num**2*nder))
ctemp.dumb
      if(logkey(ops,'m1014=debug',.false.,' ')) then
       ints=iadtwp(top)
       top=wpadti(ints+num**4)
      end if
cend.dumb
c
      call getscm(top,a,maxcor,'core')
c
c     ----- read in number of orbitals of shell, alphas and betas -----
c
      call iosys('read integer "number in shell" from rwf',nshell,
     #            a(numshl),0,' ')
      call iosys('read integer "first of shell" from rwf',nshell,
     #            a(minshl),0,' ')
      call iosys('read integer "last of shell" from rwf',nshell,
     #            a(maxshl),0,' ')
      call iosys('read real alpha from rwf',nshell**2,z(alpha),0,' ')
      call iosys('read real beta from rwf',nshell**2,z(beta),0,' ')
      call iosys('read integer "independent rotations" from rwf',nnp,
     #            a(pt),0,' ')
c
c     ----- open the integral units -----
c
      call iosys('read character "transformed integral filename"'
     $         //' from rwf',0,0,0,tints)
      call iosys('open tints as old',0,0,0,tints)
c
c     ----- and make the wa matrices -----
c
      if(logkey(ops,'m1014=debug',.false.,' ')) then
         call dumbwa(z(values),nnp,num,ntriang,z(isq),z(ds),z(wa),
     #        z(dlag),a(pt),nindep,z(alpha),z(beta),nshell,a(orbshl),
     #        a(minshl),a(maxshl),a(numshl),z(lag),z(glag),
     #        z(t1),z(t2),nder,z(ints),ops)
      else
         call makewa(z(values),nnp,num,ntriang,z(isq),z(ds),z(wa),
     #       z(dlag),a(pt),nindep,z(alpha),z(beta),nshell,a(orbshl),
     #       a(minshl),a(maxshl),a(numshl),z(lag),z(glag),
     #       z(t1),z(t2),nder,ops)
      end if
c
c     ----- and exit this link -----
c
      call iosys('close tints',0,0,0,' ')
      call chainx(0)
c
c
      stop
      end
