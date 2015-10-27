*deck @(#)kohndt.f	1.1 9/7/91
c***begin prologue     m206
c***date written       920715   (yymmdd)
c***revision date               (yymmdd)
c***keywords           m206, link 206
c***author             schneider, barry (nsf)
c***source             m206
c***purpose            massage basis function data for use
c***                   as formatted file
c***
c
c***references
c
c***routines called    iosys, util and mdutil
c***end prologue       m6010
      program m206
      implicit integer (a-z)
      character *16 bflabl(2000)
      character *4096 ops
      character *32 xform
      character *80 title
      logical logkey, prbaso, prbasn, prvec
      real *8 z
      dimension z(1)
      common a(1)
      common /io/ inp, iout
      common / memory / ioff
      equivalence (z,a)
      call drum
      call iosys ('read character options from rwf',-1,0,0,ops)
      prbaso=logkey(ops,'print=m206=old-basis',.false.,' ')
      prbasn=logkey(ops,'print=m206=new-basis',.false.,' ')
      prvec=logkey(ops,'print=m206=transformation-matrix',.false.,' ')
      write(iout,*) '                reformatting and printing basis '//
     1                               'function'
      write(iout,*) '                             and'
      write(iout,*) '                transfomation vectors from mesa'   
c----------------------------------------------------------------------c
c             get basis function information                           c
c----------------------------------------------------------------------c
      call iosys('read integer "number of basis functions" from rwf',
     $     1,nbasis,0,' ')
      if(logkey(ops,'drop',.false.,' ')) then
         idrop=1
         call iosys('read integer "old nbf" from rwf',
     $              1,oldnbf,0,' ')
         write (iout,1)
      else
         oldnbf=nbasis
         idrop=0
         write (iout,2)
      end if
c----------------------------------------------------------------------c
c                 get the basis function labels                        c
c----------------------------------------------------------------------c
      call iosys('read character "basis function labels" from rwf',
     $          -1,0,0,bflabl)
      call iosys('read integer "number of atoms" from rwf',1,nat,0,' ')
      call iosys('length of exponents on rwf',nprim,0,0,' ')
      call iosys('length of "contraction coefficients" on rwf',
     $     ncont,0,0,' ')
      call iosys('length of "number of pure functions" on rwf',
     $     ntypes,0,0,' ')
      call iosys('read integer "number basis types" from rwf',
     $     1,nbtype,0,' ')
      call iosys('length of "power of x" on rwf',lenxyz,0,0,' ')
c----------------------------------------------------------------------c
c              divide core for basis set information                   c
c----------------------------------------------------------------------c
      nnp=(nbasis+1)*nbasis/2
      nwreal=oldnbf+4*nat+ncont+nprim+nbasis*nbasis
      nwint=5*ntypes*nat+5*ntypes+3*lenxyz
      words=nwreal+iadtwp(nwint)
      call iosys ('read integer maxsiz from rwf',1,maxcor,0,' ')
      if (words.gt.maxcor) then
          call lnkerr('requested too much core')
      else
          call iosys ('write integer maxsiz to rwf',1,words,0,' ')
          call getscm(words,z,ngot,'m206',0)
      endif
c----------------------------------------------------------------------c
c            put out transformation matrix                             c
c----------------------------------------------------------------------c
      tmat=ioff
      index=tmat+nbasis*nbasis
      zan=index+oldnbf
      c=zan+nat
      cont=c+3*nat
      ex=cont+ncont
      ptprim=wpadti(ex+nprim)
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
c----------------------------------------------------------------------c
c           get basis function information                             c
c           from rwf file and put in new format                        c
c----------------------------------------------------------------------c
      call basout(z(c),z(ex),z(cont),a(ptprim),a(noprim),a(nocont),
     1            a(ptcont),nat,nprim,ntypes,nbtype,nnp,ncont,a(start),
     2            nbasis,z(zan),a(nocart),a(nobf),a(maxmom),a(mintyp),
     3            a(nx),a(ny),a(nz),a(minmom),a(index),idrop,oldnbf,
     4            ncon,prbaso,prbasn)
      call iosys('read character "transformation vector" from rwf',
     $           -1,0,0,xform)
      call iosys('read real '//xform//' from rwf',nbasis*nbasis,z(tmat),
     1            0,' ')
      if (prvec) then
          title='transformation vectors'
          call prntrm(title,z(tmat),nbasis,nbasis,nbasis,nbasis,iout)
      endif
      write(iout,*) '       transformation vector'
      ist=tmat
      do 10 i=1,nbasis
         write(iout,3) i
         write(iout,4) (z(ii),ii=ist,ist+nbasis-1)
         ist=ist+nbasis
   10 continue
      call iosys ('write integer maxsiz to rwf',1,maxcor,0,' ')
      call chainx(0)
      stop
    1 format (/,20x,'***** dropping basis functions *****')
    2 format (/,20x,'***** no basis functions dropped *****')
    3 format(/,5x,'basis function',1x,i3,/)
    4 format(5e15.8)
      end











