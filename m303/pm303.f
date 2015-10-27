*deck @(#)pm303.f	5.1  11/6/94
      subroutine pm303(z,a)
c***begin prologue     m303
c***date written       840710   (yymmdd)
c***revision date      861204   (yymmdd)
c
c         4  december 1986  pws at lanl
c                      changing 'namint' and iosys open to chracter.
c
c         17 may 1986 pws at lanl
c                      modifying to write out derivative integrals and
c                      not forming scf gradient contribution
c
c         850610       modified fairly drastically to compute
c                      effective core potential integrals (rlm,lanl).
c
c***keywords           m303, link 303, one-electron, integrals, ecp,
c                      derivatives, gradient
c***author             saxe, paul and martin, richard    (lanl)
c***source             @(#)pm303.f	5.1   11/6/94
c***purpose            computes symmetry-orbital 1-e derivative
c                      integrals over a general contraction scheme.
c***description
c     m303 currently recognizes the option strings:
c          timing         print timing statistics for this link.
c
c***references
c
c***routines called
c     m303
c       drum(mdutil)
c       iosys(io)
c       traktm(mdutil)
c       getscm(mdutil)
c       nucrep(local)
c       oneint(local)
c       chainx(mdutil)
c
c***end prologue       m303
c
      implicit integer (a-z)
      parameter (maxnbf=2000)
c
      character*4096 ops
      character*16 bflabl(maxnbf)
      character*128 rdints
      real*8 z(*)
      integer a(*)
      logical dolp
c
      parameter (nderiv=1)
c
      common /io/     inp,iout
c
      data maxcor /1/
      save maxcor
c
c     ----- collect the options string -----
c
      call iosys('read character options from rwf',-1,0,0,ops)
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
      nnp=(nbasis+1)*nbasis/2
      nbf2=nbasis*nbasis
c
c     ----- divide core for basis set information -----
c
      zan=1
      c=zan+nat
      cont=c+3*nat
      ex=cont+ncont
      s=ex+nprim
      ds=s+nnp
      if (nderiv.eq.1) then
         dhs=ds+nnp*3*nat
         ptprim=wpadti(dhs+nbf2*3*nat)
      else
         dhs=ds
         ptprim=wpadti(ds)
      end if
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
      need=nz+lenxyz
c
c     retrieve information about the most demanding shell block.
      call iosys('read integer maxprm from rwf',1,maxprm,0,' ')
      call iosys('read integer maxcont from rwf',1,mxcont,0,' ')
      call iosys('read integer maxl from rwf',1,maxl,0,' ')
      call iosys('read integer maxblk from rwf',1,maxblk,0,' ')
c
c     potential-energy integrals are most-demanding in a normal run.
c     they need:
      npint=maxprm*maxprm
      prmint=1+npint*(4*nderiv+16)+npint*(maxl+maxl+nderiv+1)
      top1=prmint+npint*(7*maxblk+9*(maxl+1)*(maxl+1)+1)
      top2=prmint+maxblk*(7*npint+mxcont*mxcont) +mxcont*maxprm
c
c     perhaps we need to do ecp's.
cps      if(dolp) then
cps         top1=max(top1,npint*((11*9*9)+7*maxblk+dlen)+310)
cps         top2=max(top2,top1+maxblk*(npint+mxcont*mxcont)+mxcont*maxprm)
cps      endif
c
      top=need+max(top1,top2)
c
      call getscm(top,z(1),ngot,'m303: main',0)
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
      call iosys('read integer "pointer to contraction '//
     $     'coefficients" from rwf',-1,a(ptcont),0,' ')
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
      call iosys('read integer "pointer to cartesians" from rwf',
     $     -1,a(mintyp),0,' ')
      call iosys('read integer "pointer to first function" from rwf',
     $     -1,a(start),0,' ')
      call iosys('read integer "power of x" from rwf',-1,a(nx),0,' ')
      call iosys('read integer "power of y" from rwf',-1,a(ny),0,' ')
      call iosys('read integer "power of z" from rwf',-1,a(nz),0,' ')
      call iosys('read character "basis function labels" from rwf',
     $          -1,0,0,bflabl)
c
c     ----- calculate one-electron integrals -----
c
      len=2*nnp*nnp
      call iosys('read character "raw derivative integral filename"'
     $         //' from rwf',0,0,0,rdints)
      call iosys('open rdints as new',len,0,0,rdints)
c
      call oneint(z(c),z(ex),z(need),z(cont),z(s),a(ptprim),a(noprim),
     #            a(nocont),a(ptcont),nat,nprim,maxcor-need+1,
     #            ntypes,nbtype,nnp,ncont,a(start),nbasis,z(zan),
     #            a(nocart),a(nobf),a(maxmom),a(mintyp),a(nx),a(ny),
     #            a(nz),a(minmom),ops,z(ds),nderiv,z(dhs))
c
c     ----- and exit gracefully -----
c
      call iosys('close rdints',0,0,0,' ')
      call chainx(0)
c
c
      stop
      end
