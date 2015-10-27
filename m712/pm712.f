*deck @(#)pm712.f	5.1 11/6/94
      subroutine pm712(z,a)
c***begin prologue     m712
c***date written       851113   (yymmdd)
c***revision date      910618   (yymmdd)
c
c   18 june    1991    rlm at lanl
c      reading group ordered density matrices from 'saoden'
c   10 january 1988    bhl at brl
c      adding capability to handle multi-reference ci gradients
c
c   05 january 1988    bhl at brl
c      bugs in expansion of group-ordered ao density
c      and in ci2pdm  removed
c
c   20 november 1987   pws at lanl
c      adding capability to read in the ao two-particle density matrix
c      directly.
c
c   11 november 1987   pws at lanl
c      adding in second derivatives.
c
c***keywords           m712, link 712, derivative integrals
c***author             saxe, paul    (lanl)
c***source             @(#)pm712.f	5.1 11/6/94
c***purpose            computes hartree-fock derivatives
c                      general contraction scheme.
c***description
c     m712 currently recognizes the option strings:
c       preexponential=n   cutoff to use on preexponential factors is
c                          10**(-n).  (default=10**-15).
c       print_grad         print the two-electron contribution the the
c                          gradient.
c       print_scf_grad     print the total scf gradient
c       timing             print timing statistics for this link.
c
c***references
c
c***routines called
c     m712
c       drum(mdutil)
c       iosys(io)
c       traktm(mdutil)
c       getscm(mdutil)
c       driver(local)
c       chainx(mdutil)
c
c***end prologue       m712
c
      implicit integer (a-z)
c
      character*4096 ops
      character*8 prtflg
      character*128 saoden
      real*8 cutexp,fpkey
      real*8 z(*)
      integer a(*)
      logical ci,mrci
      logical mcscf
      logical prnt,logkey
c
      common /io/     inp,iout
c
c
      data maxcor /1/
      data prnt/.true./
      save maxcor,prnt
c
c     ----- recover the options string -----
c
      call iosys('read character options from rwf',-1,0,0,ops)
c
c
c     ----- decide where the two-particle density matrix comes from
c
      ci=logkey(ops,'ci',.false.,' ')
c
      mcscf=.false.
      mrci=.false.
c
      if(.not.ci) then
        mcscf=logkey(ops,'mcscf',.false.,' ')
      else
        mrci=logkey(ops,'mcscf',.false.,' ')
      end if
c
c     open the density matrix unit
      if(mcscf.or.ci) then
         call iosys('read character "sao density filename" from rwf',
     $               0,0,0,saoden)
         call iosys('open saoden as old',0,0,0,saoden) 
      endif
c
c     ----- set cutoffs for preexponential and final integral -----
c
      cutexp=fpkey(ops,'gradient=preexponential',1.0d-12,' ')
c
c     see if print has been turned off externally.
      call iosys('read character "print flag" from rwf',-1,0,0,prtflg)
      if(prtflg.eq.'minimum') prnt=.false.
c
      nderiv=intkey(ops,'nderiv',1,' ')
      if(logkey(ops,'force-constants',.false.,' ')) then
         if(logkey(ops,'force-constants=numerical',.false.,' ')) then
c           numerical force constants, do nothing
         else
c           analytic force constants
            nderiv=2
         endif 
      endif
c
      if(prnt) write (iout,1) cutexp
    1 format(1x,'m712: two-electron derivative integrals',
     $     /,5x,'gradient integral pre-exponential cutoff:',1pe8.1)
c
c     ----- get lengths, dimensions, etc. needed for core allocation --
c
      call iosys('read integer "number of atoms" from rwf',1,nat,0,' ')
      call iosys('read integer "number of basis functions" from rwf',
     $     1,nbf,0,' ')
      call iosys('read integer "number of primitive functions" '//
     $     'from rwf',1,npf,0,' ')
      call iosys('length of exponents on rwf',nprim,0,0,' ')
      call iosys('length of "contraction coefficients" on rwf',
     $     ncont,0,0,' ')
      call iosys('length of "number of pure functions" on rwf',
     $     ntypes,0,0,' ')
      call iosys('read integer "number basis types" from rwf',
     $     1,nbtype,0,' ')
c..mrci
      if (mcscf) then
         ndmat=2
         nshell=3
      else if(mrci) then
         nshell=1
         ndmat=0
      else
         call iosys('read integer "number of hf density matrices" '//
     $              'from rwf',1,ndmat,0,' ')
         call iosys('read integer "number of shells" from rwf',
     $               1,nshell,0,' ')
      end if
c..mrci
      nnprim=(npf+1)*npf/2
c
c     ntypes is the total number of types used to define the lengths
c     in m301.  this total is composed of two classes:  the first
c     nbtype refer to the basis function types, the remainder refer to
c     ecp types, and are not used here.
c
      call iosys('length of "power of x" on rwf',lenxyz,0,0,' ')
      nnp=(nbf+1)*nbf/2
c
c     ----- divide core for basis set information -----
c
      if (nderiv.eq.0) then
         grad=1
         d2e=1
         zan=1
      else if (nderiv.eq.1) then
         grad=1
         d2e=1
         zan=grad+3*nat
      else if (nderiv.eq.2) then
         nd1e=3*nat
         nd2e=nd1e*(nd1e+1)/2
         grad=1
         d2e=grad+nd1e
         zan=d2e+nd2e
      end if
      c=zan+nat
      ex=c+3*nat
      cont=ex+nprim
      alpha=cont+ncont
      beta=alpha+nshell**2
      atptd=wpadti(beta+nshell**2)
      atnod=atptd+nat**2
      ptprim=atnod+nat**2
      noprim=ptprim+nat*ntypes
      ptcont=noprim+nat*ntypes
      nocont=ptcont+nat*ntypes
      mgrpno=nocont+nat*ntypes
      start=mgrpno+nat*ntypes
      pstart=start+nat*ntypes
      nocart=pstart+nat*ntypes
      nobf=nocart+ntypes
      minmom=nobf+ntypes
      maxmom=minmom+ntypes
      mintyp=maxmom+ntypes
      nx=mintyp+ntypes
      ny=nx+lenxyz
      nz=ny+lenxyz
c
      if (ci.or.mcscf) then
c
         call iosys('read integer "number of grouped 2pdm elements" '//
     $        'from rwf',1,n2pdm,0,' ')
         call iosys('read integer "number of momentum groups"  '//
     $        'from rwf',1,ngrp,0,' ')
         nnpgrp=ngrp*(ngrp+1)/2
         lentdm=30000
c
         bftgrp=nz+lenxyz
         bftcmp=bftgrp+nbf
         grpsiz=bftcmp+nbf
         gptij=grpsiz+ngrp
         gptkl=gptij+nnpgrp
         gklsiz=gptkl+nnpgrp
         tdm=iadtwp(gklsiz+nnpgrp)
         dpr=tdm+lentdm
      else
         ngrp=1
         nnpgrp=1
         n2pdm=1
         lentdm=1
         bftgrp=nz+lenxyz
         bftcmp=bftgrp
         grpsiz=bftcmp
         gptij=grpsiz
         gptkl=gptij
         gklsiz=gptkl
         tdm=iadtwp(gklsiz)
         dpr=tdm
      end if
c
      need=wpadti(dpr+nnprim*ndmat)
      rneed=iadtwp(need)
c
      if (need.gt.maxcor) then
         call getscm(need,z(1),ngot,'m712: main',0)
         maxcor=maxcor+ngot
      end if
c
c     ----- read in basis set information from the read-write file ----
c
      call iosys('read real exponents from rwf',-1,z(ex),0,' ')
      call iosys('read real "contraction coefficients" from rwf',
     $     -1,z(cont),0,' ')
      call iosys('read real coordinates from rwf',-1,z(c),0,' ')
      call iosys('read integer "pointer to primitives" from rwf',
     $     -1,a(ptprim),0,' ')
      call iosys('read integer "number of primitives" from rwf',
     $     -1,a(noprim),0,' ')
      call iosys('read integer "pointer to contraction coefficients"'//
     $     ' from rwf',-1,a(ptcont),0,' ')
      call iosys('read integer "number of contraction coefficients" '//
     $     'from rwf',-1,a(nocont),0,' ')
      call iosys('read integer "number of cartesians" from rwf',
     $     -1,a(nocart),0,' ')
      call iosys('read integer "number of pure functions" from rwf',
     $     -1,a(nobf),0,' ')
      call iosys('read integer "minimum momentum" from rwf',
     $     -1,a(minmom),0,' ')
      call iosys('read integer "maximum momentum" from rwf',
     $     -1,a(maxmom),0,' ')
      call iosys('read integer "pointer to cartesians"from rwf',
     $     -1,a(mintyp),0,' ')
      call iosys('read integer "pointer to first function" from rwf',
     $     -1,a(start),0,' ')
      call iosys('read integer "pointer to first primitive" from rwf',
     $     -1,a(pstart),0,' ')
      call iosys('read integer "power of x" from rwf',-1,a(nx),0,' ')
      call iosys('read integer "power of y" from rwf',-1,a(ny),0,' ')
      call iosys('read integer "power of z" from rwf',-1,a(nz),0,' ')
c..mrci
      if (.not.mcscf) then
         if(.not.mrci) then
            call iosys('read real alpha from rwf',
     $                  nshell**2,z(alpha),0,' ')
            call iosys('read real beta from rwf',
     $                  nshell**2,z(beta),0,' ')
         end if
      end if
c..mrci
c
      if (ci.or.mcscf) then
         call iosys('read integer "basis functions to groups" '//
     $        'from rwf',nbf,a(bftgrp),0,' ')
         call iosys('read integer "basis functions to components" '//
     $        'from rwf',nbf,a(bftcmp),0,' ')
         call iosys('read integer "momentum group size" from rwf',
     $        ngrp,a(grpsiz),0,' ')
         call iosys('read integer "group ij pointers" from rwf',
     $        nnpgrp,a(gptij),0,' ')
         call iosys('read integer "group kl pointers" from rwf',
     $        nnpgrp,a(gptkl),0,' ')
         call iosys('read integer "group pair blocksize" from rwf',
     $        nnpgrp,a(gklsiz),0,' ')
      end if
c
c     ----- calculate two electron integrals -----
c
      call driver(a(ptprim),a(noprim),nbtype,z(ex),z(c),
     #     a(nx),a(ny),a(nz),lenxyz,a(nobf),a(nocart),
     #     a(mintyp),a(minmom),a(maxmom),
     #     z(rneed),maxcor-need+1,nat,npf,nnprim,nprim,
     #     ops,cutexp,a(need),z(grad),nderiv,a(atptd),a(atnod),
     #     prnt,ndmat,z(alpha),z(beta),nshell,a(pstart),z(dpr),
     $     z(d2e),nd2e,a(bftgrp),a(bftcmp),a(grpsiz),a(gptij),
     $     a(gptkl),a(gklsiz),nbf,ngrp,nnpgrp,n2pdm,z(tdm),lentdm,
     $     ci,mcscf,mrci,z(cont),ncont,a(ptcont),a(nocont),a(start))
c
c     ----- print gradient contribution and form total gradient -----
c
      if (logkey(ops,'print=gradient=two-electron',.false.,' ')) then
         write (iout,2)
    2    format (/,5x,'the two-electron contribution to the gradients')
         call matout(z(grad),3,nat,3,nat,iout)
         if (nderiv.ge.2) then
            write (iout,3)
 3          format (/,5x,'the two-electron contribution to the force ',
     $           'constants')
            call print(z(d2e),nd2e,nd1e,iout)
         end if
      end if
c
      call iosys('read real "one-electron derivatives" from rwf',
     $     3*nat,z(zan),0,' ')
      call vadd(z(grad),z(grad),z(zan),3*nat)
c
      if (ci) then
         call iosys('write real "ci integral first derivatives" to rwf',
     $        3*nat,z(grad),0,' ')
      else if (mcscf) then
         call iosys('write real "mcscf first derivatives" to rwf',
     $              3*nat,z(grad),0,' ')
         call iosys('write real "cartesian first derivatives" to rwf',
     $              3*nat,z(grad),0,' ')
      else
         call iosys('write real "hf first derivatives" to rwf',
     $               3*nat,z(grad),0,' ')
         call iosys('write real "cartesian first derivatives" to rwf',
     $               3*nat,z(grad),0,' ')
      end if
c
      if (nderiv.ge.2) then
         call iosys('read real "one-electron force-constants" from rwf',
     $        nd2e,z(zan),0,' ')
         call vadd(z(d2e),z(d2e),z(zan),nd2e)
         call iosys('write real "integral force constants" to rwf',
     $        nd2e,z(d2e),0,' ')
      end if
c
      if (ci.and.logkey(ops,'print=gradient=integral',.false.,' ')) then
         write (iout,43)
 43      format (/5x,'ci derivatives less cphf term')
         call matout(z(grad),3,nat,3,nat,iout)
      end if
c
      if (logkey(ops,'print=gradient=total',.false.,' ')) then
         if (ci) then
         else if (mcscf) then
            write (iout,44)
 44         format (/,5x,'the mcscf first derivatives')
         else
            write (iout,4)
 4          format (/,5x,'the total scf first derivatives:')
         end if
c
         call matout(z(grad),3,nat,3,nat,iout)
         if (nderiv.ge.2) then
            write (iout,5)
 5          format (/,5x,'the force-constants less cphf contributions')
            call print(z(d2e),nd2e,nd1e,iout)
         end if
      end if
c
c
c     close the 2-particle density matrix unit.
      if(mcscf.or.ci) then
         call iosys('close saoden',0,0,0,' ')
      endif
c
c     ----- and exit gracefully -----
c
      call chainx(0)
c
c
      stop
      end
