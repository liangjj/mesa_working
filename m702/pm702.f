*deck @(#)pm702.f	5.1  11/6/94
      subroutine pm702(z,a)
c
c***begin prologue     pm702
c***date written       840710   (yymmdd)
c***revision date      871121   (yymmdd)
c
c    21 november 1987  pws at lanl
c        adding in sections for mcscf and ci gradients.
c
c    21 february 1987  pws at lanl
c        reducing memory requirements by not forming the derivative
c        integrals per se, but by directly forming the gradients.
c        also, fixing a core allocation bug per rlm...essentially
c        changing 'lenblk' to 'maxblk'.
c
c    21 august 1986        pws  at lanl
c        adding open shell capabilities. i will set it up so that the
c     code works in terms of f's, alpha's and beta's for generalized fock
c     formalism. note that the density matrices a normalized to 1.0d+00
c     and f's also go from 0 to 1.
c
c                      modified fairly drastically to compute
c                      effective core potential integrals (rlm,lanl).
c***keywords           m702, link 702, one-electron, integrals, ecp,
c                      derivatives, gradient
c***author             saxe, paul and martin, richard    (lanl)
c***source            @(#)pm702.f	5.1   11/6/94
c***purpose            computes symmetry-orbital 1-e derivative
c                      integrals over a general contraction scheme.
c***description
c     m702 currently recognizes the option strings:
c          print_s        print the overlap integrals.
c          print_t        print the kinetic energy integrals.
c          print_v        print the potential energy integrals.
c          print_lp       print the effective core potential integrals.
c          print_nuc_grad print the nuclear-repulsion contribution to the
c                         gradients.
c          print_s_grad   print the overlap contribution to the
c                         gradients.
c          print_t_grad   print the kinetic energy contribution to the
c                         gradients.
c          print_v_grad   print the potential energy contribution to the
c                         gradients.
c          print_lp_grad  print the effective core potential to the
c                         gradients.
c          print_grad     print the one-electron contribution to the
c                         gradients.
c          timing         print timing statistics for this link.
c
c***references
c
c***routines called
c     m702
c       drum(mdutil)
c       iosys(io)
c       traktm(mdutil)
c       getscm(mdutil)
c       nucrep(local)
c       oneint(local)
c       chainx(mdutil)
c
c***end prologue       pm702
c
      implicit integer (a-z)
c
      logical hf
      logical dft
      logical ci
      logical mcscf
      logical prnt
      character*4096 ops
      character*8 prtflg
      real*8 z(*)
      integer a(*)
      logical dolp,logkey
c
      common /io/     inp,iout
c
      data prnt/.true./
      data maxcor /1/
      save prnt,maxcor
c
c     ----- collect the options string -----
c
      call iosys('read character options from rwf',-1,0,0,ops)
c
      call iosys('read character "print flag" from rwf',-1,0,0,prtflg)
      if(prtflg.eq.'minimum') prnt=.false.
c
      if(prnt) then
         write(iout,10001)
10001    format(1x,'m702: one-electron derivative integrals ')
      endif
c
c     ----- decide on type -----
c
      hf=logkey(ops,'hf',.false.,' ')
      dft=logkey(ops,'dft',.false.,' ')
c     --- for now, signal dft by making hf .true.
      if(dft) hf=.true.
      ci=logkey(ops,'ci',.false.,' ')
c
      mcscf=.false.
      if(.not.ci) then
       mcscf=logkey(ops,'mcscf',.false.,' ')
      end if
      if((.not.hf).and.(.not.dft).and.(.not.mcscf).and.(.not.ci)) then
         call lnkerr('cannot determine run type.'
     $              //' specify hf,dft,mcscf,ci in route')
      endif
c
c     ----- get the lengths of arrays needed for core allocation -----
c     ntypes is the total number of types used to define the lengths
c     in m102.  this total is composed of two sets: the first nbtype
c     types are used to mark the basis function types.  the remainder
c     refer to the ecp types.
c
      call iosys('read integer "number of atoms" from rwf',1,nat,0,' ')
      call iosys('read integer "number of basis functions" from rwf',
     $     1,nbasis,0,' ')
      call iosys('length of exponents on rwf',nprim,0,0,' ')
      call iosys('length of "contraction coefficients" on rwf',
     $     ncont,0,0,' ')
      call iosys('length of "number of pure functions" on rwf',
     $     ntypes,0,0,' ')
      call iosys('read integer "number basis types" from rwf',
     $     1,nbtype,0,' ')
      call iosys('length of "power of x" on rwf',lenxyz,0,0,' ')
      if (mcscf.or.ci) then
         ndmat=1
      else
         call iosys('read integer "number of hf density matrices" '//
     $        'from rwf',1,ndmat,0,' ')
      end if
      nnp=(nbasis+1)*nbasis/2
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
      nd1e=3*nat
      nd2e=nd1e*(nd1e+1)/2
c
c----------------------------------------------------------------------
c ndmat is the number of density matrices, which is one less than the
c       number of orbital types. (virtual orbitals do not contribute)
c----------------------------------------------------------------------
c
c     ----- divide core for basis set information -----
c
      f=1
      d=f+ndmat+1
      d2e=d+nnp
      if (nderiv.eq.1) then
         ld2e=d2e
         grad=d2e
      else if (nderiv.eq.2) then
         ld2e=d2e+nd2e
         grad=ld2e+nd2e
      end if
      tempg=grad+3*nat
      zan=tempg+3*nat
      c=zan+nat
      cont=c+3*nat
      ex=cont+ncont
      s=ex+nprim
      ptprim=wpadti(s+nnp)
      noprim=ptprim+ntypes*nat
      nocont=noprim+ntypes*nat
      ptcont=nocont+ntypes*nat
      start=ptcont+ntypes*nat
      nocart=start+ntypes*nat
      nobf=nocart+ntypes
      minmom=nobf+ntypes
      maxmom=minmom+ntypes
      mintyp=maxmom+ntypes
      nx=mintyp+ntypes
      ny=nx+lenxyz
      nz=ny+lenxyz
      need=iadtwp(nz+lenxyz)
c
c     retrieve information about the most demanding shell block.
c
      call iosys('read integer maxprm from rwf',1,maxprm,0,' ')
      call iosys('read integer maxcont from rwf',1,mxcont,0,' ')
      call iosys('read integer maxl from rwf',1,maxl,0,' ')
      call iosys('read integer maxblk from rwf',1,maxblk,0,' ')
      call iosys('read integer dolp from rwf',1,dolp,0,' ')
      call iosys('read integer d1maxblk from rwf',1,dlen,0,' ')
      call iosys('read integer d2maxblk from rwf',1,d2len,0,' ')
      if(nderiv.eq.2) dlen=d2len
c
c     potential-energy integrals are most-demanding in a normal run.
c     they need:
c
      npint=maxprm*maxprm
      prmint=1+npint*(4*nderiv+16)+npint*(maxl+maxl+nderiv+1)
      top1=prmint+npint*(7*maxblk+9*(maxl+1)*(maxl+1)+1)
      top2=prmint+maxblk*(7*npint+mxcont*mxcont) +mxcont*maxprm
c
c     perhaps we need to do ecp's.
c
      if(dolp) then
         if(nderiv.lt.2) then
            top1=max(top1,npint*((11*9*9)+7*maxblk+dlen)+310)
         else if(nderiv.eq.2) then
            top1=max(top1,npint*((11*9*9)+28*maxblk+dlen)+310)
         endif
         top2=max(top2,top1+maxblk*(npint+mxcont*mxcont)+mxcont*maxprm)
      endif
      top=need+max(top1,top2)
c
      if (top.gt.maxcor) then
         call getscm(wpadti(top),z(1),ngot,'m702: main',0)
         maxcor=iadtwp(ngot)
      end if
c
c     ----- read in basis set information from read-write file -----
c
      call iosys('read real exponents from rwf',-1,z(ex),0,' ')
      call iosys('read real "contraction coefficients" from rwf',
     $     -1,z(cont),0,' ')
      call iosys('read real "nuclear charges" from rwf',
     $     -1,z(zan),0,' ')
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
      call iosys('read integer "power of x" from rwf',-1,a(nx),0,' ')
      call iosys('read integer "power of y" from rwf',-1,a(ny),0,' ')
      call iosys('read integer "power of z" from rwf',-1,a(nz),0,' ')
      if (hf) call iosys('read real f from rwf',ndmat+1,z(f),0,' ')
c
c     ----- calculate the nuclear repulsion energy -----
c
      call nucrep(z(zan),z(c),nat,z(grad))
      if (nderiv.ge.2) call nuctwo(z(zan),z(c),nat,z(d2e),nd2e)
c
      if (logkey(ops,'print=gradient=nuclear-repulsion',.false.,' '))
     $     then
         write (iout,1)
    1    format ('1',//,'     the nuclear-repulsion contribution to',
     #          //' the gradients',//)
         call matout(z(grad),3,nat,3,nat,iout)
c
         if (nderiv.eq.2) then
            write (iout,201)
  201       format (//'      the nuclear-repulsion contribution to ',
     #              'the force constants')
            call print(z(d2e),nd2e,nd1e,iout)
         end if
      end if
c
      if(logkey(ops,'nacme',.false.,' ')) then
         write(iout,1000)
 1000    format(5x,'nacme run')
         call rzero(z(grad),3*nat)
      end if
c
c     ----- calculate one-electron integrals -----
c
      call oneint(z(c),z(ex),z(need),z(cont),z(s),a(ptprim),a(noprim),
     #            a(nocont),a(ptcont),nat,nprim,maxcor-need+1,
     #            ntypes,nbtype,nnp,ncont,a(start),nbasis,z(zan),
     #            a(nocart),a(nobf),a(maxmom),a(mintyp),a(nx),a(ny),
     #            a(nz),a(minmom),ops,nderiv,z(grad),z(tempg),
     #            z(d),z(f),ndmat,z(d2e),z(ld2e),nd2e,nd1e,hf,ci,mcscf,
     $     dft)
c
c     ----- save the gradients to the rwf -----
c
      call iosys('write real "one-electron derivatives" to rwf',
     $     3*nat,z(grad),0,' ')
      if (nderiv.ge.2) then
         call iosys('write real "one-electron force-constants" to rwf',
     $        nd2e,z(d2e),0,' ')
      end if
c
      if (logkey(ops,'print=gradient=one-electron',.false.,' ')) then
         write (iout,20)
 20      format ('1',//,'     the nuclear-repulsion, overlap and ',
     $        'one-electron contribution to the gradients')
         call matout(z(grad),3,nat,3,nat,iout)
c
         if (nderiv.ge.2) then
            write (iout,21)
 21         format (//'      the nuclear-repulsion, overlap and ',
     $           'one-electron contribution to the force constants')
            call print(z(d2e),nd2e,nd1e,iout)
         end if
      end if
c
c     ----- print statistics if desired -----
c
c
c     ----- and exit gracefully -----
c
      call chainx(0)
c
c
      return
      end
