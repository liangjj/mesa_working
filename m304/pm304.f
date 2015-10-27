*deck @(#)pm304.f	5.1  11/6/94
      subroutine pm304(z,a)
c***begin prologue     m304
c***date written       911204   (yymmdd)
c***revision date         (yymmdd)
c                      
c                     
c***keywords           m304, link 304, one-electron, spin-orbit
c***author             m. braunstein (lanl) (modified from m302)
c***source             
c***purpose            computes 1-e spin-orbit
c                      integrals over a
c                      general contraction scheme.
c***description
c     m304 currently recognizes the option strings:
c       print_so    print the spin-orbit integrals
c       timing      print timing statistics for this link.
c
c***references
c
c***routines called
c     m304
c       drum(mdutil)
c       iosys(io)
c       traktm(mdutil)
c       getcon(io)
c       getscm(mdutil)
c       soint(local)
c       chainx(mdutil)
c
c***end prologue       m304
      implicit integer (a-z)
      parameter (maxnbf=2000)
c
      character*4096 ops
      character*128 namint
      character*16 bflabl(maxnbf)
      logical dolp,logkey
      real*8 z(*)
      integer a(*)
c
      common /io/     inp,iout
c
      data maxcor /1/
      data bufsiz/4096/
c
  100 format(1x,'m304: skip one-electron spin-orbit integrals')
c     ----- collect the options string -----
c
      call iosys('read character options from rwf',-1,0,0,ops)
c
c     see if we are to compute integrals.
      if(logkey(ops,'int=reuseso',.false.,' ').or.
     $   logkey(ops,'noints',.false.,' ')) then
         write(iout,100)
c        restore the rwf file.
c        get some buffer space.
         call getscm(bufsiz,z(1),ngot,'m304: buffer',0)
         call iosys('read character "integral filename" from rwf',
     $               0,0,0,namint)
         call iosys('open ints as old',0,0,0,namint)
         call iosys('copy "asox integrals" from ints to rwf',
     $               bufsiz,a,0,'noerror')
         call iosys('copy "asoy integrals" from ints to rwf',
     $               bufsiz,a,0,'noerror')
         call iosys('copy "asoz integrals" from ints to rwf',
     $               bufsiz,a,0,'noerror')
      else
c
c        ----- get the lengths of arrays needed for core allocation -----
c        ntypes is the total number of types used to define the lengths
c        in m102.  this total is composed of two sets: the first nbtype
c        types are used to mark the basis function types.  the remainder
c        refer to the ecp types.
c
         call iosys('read integer "number of atoms" from rwf',1,nat,
     $               0,' ')
         call iosys('read integer "number of basis functions" from rwf',
     $                  1,nbasis,0,' ')
         call iosys('length of exponents on rwf',nprim,0,0,' ')
         call iosys('length of "contraction coefficients" on rwf',
     $               ncont,0,0,' ')
         call iosys('length of "number of pure functions" on rwf',
     $               ntypes,0,0,' ')
         call iosys('read integer "number basis types" from rwf',
     $               1,nbtype,0,' ')
         call iosys('length of "power of x" on rwf',lenxyz,0,0,' ')
         nnp=(nbasis+1)*nbasis/2
c
c        ----- divide core for basis set information -----
c
         zan=1
         c=zan+nat
         cont=c+3*nat
         ex=cont+ncont
         s=ex+nprim
         zeff=s+nbasis*nbasis
         ptprim=wpadti(zeff+nat*ntypes)
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
         rneed=iadtwp(need)
c
c        allocate core for integral computation.
c        retrieve information about the most demanding shell block.
         call iosys('read integer maxprm from rwf',1,maxprm,0,' ')
         call iosys('read integer maxcont from rwf',1,mxcont,0,' ')
         call iosys('read integer maxl from rwf',1,maxl,0,' ')
         call iosys('read integer maxblk from rwf',1,maxblk,0,' ')
         call iosys('read integer dolp from rwf',1,dolp,0,' ')
c
c        potential-energy integrals are most-demanding in a normal run.
c        they need:
         npint=maxprm*maxprm
         prmint=1+16*npint+2*npint*(maxl+1)
         top1=prmint+npint*(maxblk+6*(maxl+1)*(maxl+1)+1)
         top1=max(top1,prmint+nbasis*nbasis)
         top2=prmint+maxblk*(npint+mxcont*mxcont) +mxcont*maxprm
c        effective core potentials?
         if(dolp) then
            top1=max(top1,npint*(343+maxblk)+310)
         endif
         top=need+wpadti(max(top1,top2))
c
         if (top.gt.maxcor) then
            call getscm(top,z(1),maxcor,'m304: main',0)
            maxcor=iadtwp(maxcor)
         end if
c
c        ----- read in basis set information from read-write file -----
c
         call iosys('read real exponents from rwf',-1,z(ex),0,' ')
         call iosys('read real "contraction coefficients" from rwf',
     $               -1,z(cont),0,' ')
         call iosys('read real "nuclear charges" from rwf',
     $               -1,z(zan),0,' ')
         call iosys('read real coordinates from rwf',-1,z(c),0,' ')
         call iosys('read integer "pointer to primitives" from rwf',
     $              -1,a(ptprim),0,' ')
         call iosys('read integer "number of primitives" from rwf',
     $              -1,a(noprim),0,' ')
         call iosys('read integer "pointer to contraction '//
     $              'coefficients" from rwf',-1,a(ptcont),0,' ')
         call iosys('read integer "number of contraction coefficients"'
     $              //' from rwf',-1,a(nocont),0,' ')
         call iosys('read integer "number of cartesians" from rwf',
     $              -1,a(nocart),0,' ')
         call iosys('read integer "number of pure functions" from rwf',
     $              -1,a(nobf),0,' ')
         call iosys('read integer "minimum momentum" from rwf',
     $               -1,a(minmom),0,' ')
         call iosys('read integer "maximum momentum" from rwf',
     $              -1,a(maxmom),0,' ')
         call iosys('read integer "pointer to cartesians" from rwf',
     $               -1,a(mintyp),0,' ')
         call iosys('read integer "pointer to first function" from rwf',
     $               -1,a(start),0,' ')
         call iosys('read integer "power of x" from rwf',-1,a(nx),0,' ')
         call iosys('read integer "power of y" from rwf',-1,a(ny),0,' ')
         call iosys('read integer "power of z" from rwf',-1,a(nz),0,' ')
         call iosys('read character "basis function labels" from rwf',
     $               -1,0,0,bflabl)
c
c        ----- calculate the atomic spin-orbit integrals -----
c
         call soint(z(c),z(ex),z(rneed),a(need),z(cont),z(s),a(ptprim),
     #               a(noprim),a(nocont),a(ptcont),nat,nprim,
     #               maxcor-rneed+1,ntypes,nbtype,nnp,ncont,a(start),
     #               nbasis,z(zan),a(nocart),a(nobf),a(maxmom),
     #               a(mintyp),a(nx),a(ny),a(nz),a(minmom),dolp,
     #               bflabl,ops,z(zeff))
c
      endif
c     ----- and exit gracefully -----
c
      call chainx(0)
c
c
      stop
      end
