*deck @(#)pm323.f	5.1  11/6/94
      subroutine pm323(z,a)
c***begin prologue     m323
c***date written       891208   (yymmdd)
c
c***keywords           m323, link 323, derivative integrals
c***author             saxe, paul (lanl) and lengsfield, byron (llnl)
c***source             @(#)pm323.f	5.1   11/6/94
c***purpose            computes derivative integrals for a
c                      general contraction scheme and writes
c                      chained integral file.
c***description
c     m323 is a modified version of m313
c     m323 currently recognizes the option strings:
c       preexponential=n   cutoff to use on preexponential factors is
c                          10**(-n).  (default=10**-15).
c       cutoff=n           cutoff for retaining derivative integrals
c                          10**(-n).  (default=10**-10).
c
c***references
c
c***routines called
c       drum(mdutil)
c       iosys(io)
c       traktm(mdutil)
c       getscm(mdutil)
c       driver(local)
c       chainx(mdutil)
c
c***end prologue       m323
c
      implicit integer (a-z)
c
      character*4096 ops
      character*8 prtflg
      character*128 rdints
      real*8 cutoff,cutexp
      integer a(*)
      real*8 z(*)
      logical keyval,prnt
c
      common /io/     inp,iout
      common /toler/  cutoff
c
c
      data maxcor /1/
      data prnt/.true./
      parameter (nderiv=1)
      save maxcor,prnt
c
c     ----- recover the options string -----
c
      call iosys('read character options from rwf',-1,0,0,ops)
c
c     ----- set cutoffs for preexponential and final integral -----
c
      if (.not.keyval(ops,'m323=preex',n)) n=13
      cutexp=10.0d+00**(-n)
      if (.not.keyval(ops,'m323=cutoff',n)) n=9
      cutoff=10.0d+00**(-n)
c
c     see if print has been turned off externally.
      call iosys('read character "print flag" from rwf',-1,0,0,prtflg)
      if(prtflg.eq.'minimum') prnt=.false.
c
      maxbuf=intkey(ops,'m323=maxbuf',500000,' ')
c
      minbuf=intkey(ops,'m323=minbuf',10000,' ')
c
      call iosys('read integer "number of atoms" from rwf',1,nat,0,' ')
c
      nder=3*nat
      lenbuf=maxbuf/(2*nder)
      lenbuf=min(lenbuf,50000)
      if(lenbuf.lt.minbuf) then
         write(iout,111) maxbuf,lenbuf,minbuf,3*nat
 111     format(/,'  minimum buffer size failure ',/,
     $'              maxbuf lenbuf minbuf nat3 ',4i8)
         call lnkerr(' m323: buffer size error ')
      endif
c
      if(prnt) then
         write (iout,1) cutexp,cutoff
    1    format(1x,'m323:',
     $        /,5x,'gradient integral pre-exponential cutoff:',1pe8.1,
     $        /,5x,'gradient integral retention       cutoff:',1pe8.1)
      end if
c
c     ----- get lengths, dimensions, etc. needed for core allocation --
c
      call iosys('read integer "number of basis functions" from rwf',
     $     1,nbf,0,' ')
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
c
      nnp=(nbf+1)*nbf/2
c
c     ----- open the integrals file -----
c
      call iosys('read character "raw derivative integral filename"'
     $     //' from rwf',0,0,0,rdints)
      call iosys('open rdints as old',0,0,0,rdints)
c
      call iosys('write integer "derivative buffer length" to rdints',
     $     1,lenbuf,0,' ')
      call iosys('write integer "number of atoms" to rdints',
     $           1,nat,0,' ')
      call iosys('create real file "unsorted derivative integrals" '//
     $     'on rdints',-1,0,0,' ')
c
c     ----- divide core for basis set information -----
c
      c=1
      ex=c+3*nat
      cont=ex+nprim
      ints=cont+ncont
      ptprim=wpadti(ints+lenbuf*3*nat)
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
      numbuf=nz+lenxyz
      last=numbuf+3*nat
      labels=last+3*nat
c
c     ---- on 32 bit integer / 64 bit real machines, 'labels'
c            must be a 2 x lenbuf array
c
      if (wptoin(1).eq.2) then
         need=labels+2*lenbuf*3*nat
         lenb2=2*lenbuf
      else
         need=labels+lenbuf*3*nat
         lenb2=lenbuf
      end if
      rneed=iadtwp(need)
c
c     ---- ensure that need and rneed point to the same location -----
c
      need=wpadti(rneed)
c
      call getscm(need,z(1),ngot,'m323: main',0)
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
     #            z(ints),a(labels),lenbuf,nbf,nnp,cutoff,
     #            a(last),a(numbuf),lenb2)
c
c     ----- and exit gracefully -----
c
      call iosys('close rdints',0,0,0,' ')
      call chainx(0)
c
c
      stop
      end
