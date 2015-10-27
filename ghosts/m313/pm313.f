*deck @(#)pm313.f	1.1  11/30/90
      subroutine pm313(z,a)
c
c***begin prologue     m313
c***date written       851113   (yymmdd)
c***revision date      861204   (yymmdd)
c
c   4 december 1986  pws at lanl
c        changing 'namint' and iosys open to character.
c
c***keywords           m313, link 313, derivative integrals
c***author             saxe, paul    (lanl)
c***source             @(#)pm313.f	1.1   11/30/90
c***purpose            computes derivative integrals for a
c                      general contraction scheme.
c***description
c     m313 currently recognizes the option strings:
c       preexponential=n   cutoff to use on preexponential factors is
c                          10**(-n).  (default=10**-15).
c
c***references
c
c***routines called
c     m313
c       drum(mdutil)
c       iosys(io)
c       traktm(mdutil)
c       getscm(mdutil)
c       driver(local)
c       chainx(mdutil)
c
c***end prologue       m313
c
      implicit integer (a-z)
c
      character*4096 ops
      character*8 scratch,prtflg
      character*128 namint
      real*8 cutoff,cutexp
      integer a(*)
      real*8 z(*)
      logical keyval,prnt
c
c..bhl..unicos      common //       z(1)
      common /io/     inp,iout
      common /toler/  cutoff
c
c..bhl..unicos      equivalence (a,z)
c
      data maxcor /1/
      data lenbuf /10000/
      data prnt/.true./
c
c     ----- open files and retrieve important data -----
c
c..bhl..unicos      call drum
c
c     ----- recover the options string -----
c
      call iosys('read character options from rwf',-1,0,0,ops)
c
c     ----- start timing routines if enabled -----
c
c
c     ----- set cutoffs for preexponential and final integral -----
c
      if (.not.keyval(ops,'preex',n)) n=15
      cutexp=10.0d+00**(-n)
      if (.not.keyval(ops,'cutoff',n)) n=10
      cutoff=10.0d+00**(-n)
c
c     see if print has been turned off externally.
      call iosys('read character "print flag" from rwf',-1,0,0,prtflg)
      if(prtflg.eq.'minimum') prnt=.false.
c
      nderiv=1
c
      if(prnt) write (iout,1) cutexp
    1 format(1x,'m313:',
     $     /,5x,'gradient integral pre-exponential cutoff:',1pe8.1)
c
c     ----- get lengths, dimensions, etc. needed for core allocation --
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
c
c     ntypes is the total number of types used to define the lengths
c     in m301.  this total is composed of two classes:  the first
c     nbtype refer to the basis function types, the remainder refer to
c     ecp types, and are not used here.
c
      call iosys('length of "power of x" on rwf',lenxyz,0,0,' ')
      nnp=(nbf+1)*nbf/2
c     nprim=intowp(nprim)
c     ncont=intowp(ncont)
c
c     ----- open the integrals file -----
c
      call iosys('read character "integral filename" from rwf',
     $     0,0,0,namint)
      call iosys('open ints as old',0,0,0,namint)
      call iosys('write integer "derivative buffer length" to ints',
     $     1,lenbuf,0,' ')
      call iosys('write integer "number of atoms" to ints',1,nat,0,' ')
      call iosys('create real file "unsorted derivative integrals" '//
     $     'on ints',-1,0,0,' ')
c
c     ----- divide core for basis set information -----
c
      c=1
      ex=c+3*nat
      cont=ex+nprim
      ints=cont+ncont
      ptprim=wpadti(ints+lenbuf)
      noprim=ptprim+nat*ntypes
      ptcont=noprim+nat*ntypes
      nocont=ptcont+nat*ntypes
      start=nocont+nat*ntypes
      nocart=start+nat*ntypes
      nobf=nocart+ntypes
      minmom=nobf+ntypes
      maxmom=minmom+ntypes
      mintyp=maxmom+ntypes
      nx=mintyp+ntypes
      ny=nx+lenxyz
      nz=ny+lenxyz
      labels=nz+lenxyz
c
c     ---- on 32 bit integer / 64 bit real machines, 'labels'
c            must be a 2 x lenbuf array
c
      if (wptoin(1).eq.2) then
         need=labels+2*lenbuf
      else
         need=labels+lenbuf
      end if
      rneed=iadtwp(need)
c
c     ---- ensure that need and rneed point to the same location -----
c
      need=wpadti(rneed)
c
      call getscm(need,z(1),ngot,'m313: main',0)
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
      call iosys('read integer "power of x" from rwf',-1,a(nx),0,' ')
      call iosys('read integer "power of y" from rwf',-1,a(ny),0,' ')
      call iosys('read integer "power of z" from rwf',-1,a(nz),0,' ')
c
c     ----- calculate two electron integrals -----
c
      call driver(a(ptprim),a(noprim),nbtype,z(ex),
     #            z(c),
     #            a(nx),a(ny),a(nz),lenxyz,a(nobf),a(nocart),
     #            a(mintyp),a(minmom),a(maxmom),
     #            z(rneed),nat,nprim,nderiv,
     #            ops,cutexp,a(need),
     #            prnt,a(nocont),a(ptcont),z(cont),ncont,a(start),
     #            z(ints),a(labels),lenbuf,nbf,nnp,cutoff)
c
c     ----- print statistics if desired -----
c
c
c     ----- and exit gracefully -----
c
      call chainx(0)
c
c
      stop
      end
