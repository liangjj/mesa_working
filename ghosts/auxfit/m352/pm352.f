*deck %W%  %G%
      subroutine pm352(z,a)
c***begin prologue     m352
c***date written       840710   (yymmdd)
c***revision date      920417   (yymmdd)
c   april 17,    1992  rlm at lanl
c                      modified fairly drastically to compute
c                      3-center overlap and two-electron integrals 
c                      for auxiliary basis.
c***keywords           m352, link 352, one-electron, integrals, ecp
c***author             saxe, paul and martin, richard    (lanl)
c***source             %W%   %G%
c***purpose            computes symmetry-orbital 1-e integrals over a
c                      general contraction scheme.
c***description
c
c***references
c
c***routines called
c     pm352
c       drum(mdutil)
c       iosys(io)
c       getscm(mdutil)
c       nucrep(local)
c
c***end prologue       m352
      implicit integer (a-z)
      parameter (maxnbf=2000)
c
      character*4096 ops
      character*128 namint
      character*16 bflabl(maxnbf)
      character*16 xbflbl(maxnbf)
      logical logkey,debug
      real*8 z(*)
      integer a(*)
c
      parameter (debug=.false.)
c
      common /io/     inp,iout
c
      data maxcor /1/
      data bufsiz/4096/
      save maxcor,bufsiz
c
  100 format(1x,'m352: skip one-electron integrals')
  105 format(5x,'auxiliary 3-center overlap matrix:')
  110 format(5x,'auxiliary basis overlap matrix')
  120 format(5x,'auxiliary two-electron integrals')
  130 format(5x,'projected auxiliary two-electron integrals:')
c
c     ----- collect the options string -----
c
      call iosys('read character options from rwf',-1,0,0,ops)
c
c     see if we are to compute integrals.
      if(logkey(ops,'int=reuse',.false.,' ').or.
     $   logkey(ops,'int=reuse1',.false.,' ').or.
     $   logkey(ops,'noints',.false.,' ')) then
         write(iout,100)
c        restore the rwf file.
c        get some buffer space.
         call getscm(bufsiz,z(1),ngot,'m352: buffer',0)
         call iosys('read character "integral filename" from rwf',
     $               0,0,0,namint)
         call iosys('open ints as old',0,0,0,namint)
         call iosys('copy "nuclear repulsion energy" from ints to rwf',
     $               bufsiz,a,0,'noerror')
         call iosys('copy "overlap integrals" from ints to rwf',
     $               bufsiz,a,0,'noerror')
         call iosys('copy "kinetic integrals" from ints to rwf',
     $               bufsiz,a,0,'noerror')
         call iosys('copy "potential integrals" from ints to rwf',
     $               bufsiz,a,0,'noerror')
      else
c
c        ----- get the lengths of arrays needed for core allocation -----
c        ntypes is the total number of types used to define the lengths
c        in m102.  this total is composed of two sets: the first nbtype
c        types are used to mark the basis function types.  the remainder
c        refer to the ecp types.
c
         call iosys('read integer "number of atoms" from rwf',
     $               1,nat,0,' ')
         call iosys('read integer "number of basis functions" from rwf',
     $               1,nbasis,0,' ')
         call iosys('read integer "number of auxiliary basis functions"'
     $              //' from rwf',1,nbasx,0,' ')
         call iosys('length of exponents on rwf',nprim,0,0,' ')
         call iosys('length of "auxiliary exponents" on rwf',
     $               npaux,0,0,' ')
         call iosys('length of "contraction coefficients" on rwf',
     $               ncont,0,0,' ')
         call iosys('length of "auxiliary contraction coefficients"'
     $               //' on rwf',ncontx,0,0,' ')
         call iosys('length of "number of pure functions" on rwf',
     $               ntypes,0,0,' ')
         call iosys('read integer "number basis types" from rwf',
     $               1,nbtype,0,' ')
         call iosys('length of "power of x" on rwf',lenxyz,0,0,' ')
         nnp=(nbasis+1)*nbasis/2
         nnpx=(nbasx+1)*nbasx/2
c
c        ----- divide core for basis set information -----
c
         zan=1
         c=zan+nat
         cont=c+3*nat
         xcont=cont+ncont
         ex=xcont+ncontx
         exaux=ex+nprim
         s=exaux+npaux
         sx=s+nnp*nbasx
         sinvx=sx+nnpx
         u=sinvx+nnpx
         t1=u+nbasx*nbasx
         t2=t1+nbasx*nbasx
         eigval=t2+nbasx*nbasx
         triang=eigval+nbasx
         ptprim=wpadti(triang+nnpx)
         xptprm=ptprim+ntypes*nat
         noprim=xptprm+ntypes*nat
         xnoprm=noprim+ntypes*nat
         nocont=xnoprm+ntypes*nat
         xnocon=nocont+ntypes*nat
         ptcont=xnocon+ntypes*nat
         xptcon=ptcont+ntypes*nat
         start=xptcon+ntypes*nat
         xstart=start+ntypes*nat
         nocart=xstart+ntypes*nat
         nobf=nocart+ntypes
         minmom=nobf+ntypes
         maxmom=minmom+ntypes
         mintyp=maxmom+ntypes
         nx=mintyp+ntypes
         ny=nx+lenxyz
         nz=ny+lenxyz
c
c        allocate core for integral computation.
c        retrieve information about the most demanding shell block.
         call iosys('read integer maxprm from rwf',1,maxprm,0,' ')
         call iosys('read integer "auxiliary maxprm" from rwf',
     $               1,amxprm,0,' ')
         call iosys('read integer maxcont from rwf',1,mxcont,0,' ')
         call iosys('read integer "auxiliary maxcont" from rwf',
     $               1,amxcon,0,' ')
         call iosys('read integer maxl from rwf',1,maxl,0,' ')
         call iosys('read integer "auxiliary maxl" from rwf',
     $               1,amaxl,0,' ')
         call iosys('read integer maxblk from rwf',1,maxblk,0,' ')
         call iosys('read integer "auxiliary maxblk" from rwf',
     $               1,amxblk,0,' ')
c
c some convenient quantities
c
         mmax=2*amaxl
         mindx=((mmax+1)*(mmax+2))/2

c
c allocate some more stuff
c
         tmpsiz=nz+lenxyz
         fna0=tmpsiz+mmax+1
         lval=fna0+mmax+1
         mycart=lval+(mmax+1)*mindx*3
         yn0=mycart+mmax+1
         need=yn0+mmax+1
         rneed=iadtwp(need)

c
c        we need:
         npint=maxprm*maxprm*amxprm
         top1=1+9*npint
         top2=1+npint*maxblk*amxblk
         top=max(top1,top2)
         top1=top+3*npint*(maxl+1)*(maxl+1)*(amaxl+1)
         top2=1+maxblk*amxblk*(mxcont*mxcont*amxcon)
         top=max(top1,top2)+mxcont*amxprm
         top=need+wpadti(top)
c
         if (top.gt.maxcor) then
            call getscm(top,z(1),maxcor,'m352: main',0)
            maxcor=iadtwp(maxcor)
         end if
c
c        ----- read in basis set information from read-write file -----
c
         call iosys('read real exponents from rwf',-1,z(ex),0,' ')
         call iosys('read real "auxiliary exponents" from rwf',
     $               -1,z(exaux),0,' ')
         call iosys('read real "contraction coefficients" from rwf',
     $               -1,z(cont),0,' ')
         call iosys('read real "auxiliary contraction coefficients"'
     $               //' from rwf',-1,z(xcont),0,' ')
         call iosys('read real "nuclear charges" from rwf',
     $               -1,z(zan),0,' ')
         call iosys('read real coordinates from rwf',-1,z(c),0,' ')
         call iosys('read integer "pointer to primitives" from rwf',
     $              -1,a(ptprim),0,' ')
         call iosys('read integer "pointer to auxiliary primitives"'
     $              //' from rwf',-1,a(xptprm),0,' ')
         call iosys('read integer "number of auxiliary primitives"'
     $              //' from rwf',-1,a(xnoprm),0,' ')
         call iosys('read integer "number of primitives" from rwf',
     $              -1,a(noprim),0,' ')
         call iosys('read integer "pointer to contraction '//
     $              'coefficients" from rwf',-1,a(ptcont),0,' ')
         call iosys('read integer "pointer to auxiliary contraction '//
     $              'coefficients" from rwf',-1,a(xptcon),0,' ')
         call iosys('read integer "number of contraction coefficients"'
     $              //' from rwf',-1,a(nocont),0,' ')
         call iosys('read integer "number of auxiliary contraction'
     $              //' coefficients" from rwf',-1,a(xnocon),0,' ')
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
         call iosys('read integer "pointer to first auxiliary function"'
     $               //' from rwf',-1,a(xstart),0,' ')
         call iosys('read integer "power of x" from rwf',-1,a(nx),0,' ')
         call iosys('read integer "power of y" from rwf',-1,a(ny),0,' ')
         call iosys('read integer "power of z" from rwf',-1,a(nz),0,' ')
         call iosys('read character "basis function labels" from rwf',
     $               -1,0,0,bflabl)
         call iosys('read character "auxiliary basis function labels"'
     $               //' from rwf',-1,0,0,xbflbl)
c
c        ----- calculate three-center overlap integrals -----
c
         call centr3(z(c),z(ex),z(rneed),a(need),z(cont),z(s),a(ptprim),
     #               a(noprim),a(nocont),a(ptcont),nat,nprim,
     #               maxcor-rneed+1,ntypes,nbtype,nnp,ncont,a(start),
     #               nbasis,z(zan),a(nocart),a(nobf),a(maxmom),
     #               a(mintyp),a(nx),a(ny),a(nz),a(minmom),
     #               ops,z(exaux),npaux,a(xptprm),a(xnoprm),
     $               a(xnocon),a(xptcon),a(xstart),z(xcont),nbasx,
     $               ncontx)
c
c        ----- perhaps print the integrals and store -----
         if(logkey(ops,'print=int=s3aux',.false.,' ')) then
            write(iout,105)
            pts=0
            do 200 k=1,nbasx
               write(iout,*) 'auxiliary function k:',k,xbflbl(k)
               call trtosq(z(t1),z(s+pts),nbasis,nnp)
               call wlmat(z(t1),nbasis,nbasis,bflabl,bflabl)
c               pts=pts+nbasx
               pts=pts+nnp
  200       continue
         end if
         call iosys('write real "expansion overlap integrals" on rwf',
     $               nnp*nbasx,z(s),0,' ')
c
c        ----- calculate the auxiliary basis overlap matrix
         call overlp(z(c),z(exaux),z(rneed),a(need),z(xcont),z(sx),
     $               a(xptprm),a(xnoprm),a(xnocon),a(xptcon),nat,npaux,
     $               maxcor-rneed+1,ntypes,nbtype,nnpx,ncontx,a(xstart),
     $               nbasx,a(nocart),a(nobf),a(maxmom),a(mintyp),
     $               a(nx),a(ny),a(nz),a(minmom),ops)
         if(logkey(ops,'print=int=auxs',.false.,' ')) then
            write(iout,110)
            call trtosq(z(t1),z(sx),nbasx,nnpx)
            call wlmat(z(t1),nbasx,nbasx,xbflbl,xbflbl)
         end if
c
c        ----- calculate the inverse auxiliary overlap matrix -----
         iprint=0
         call sinvrt(z(sx),z(sinvx),z(u),z(eigval),z(t1),z(t2),
     $               nbasx,nnpx,z(triang),iprint)
c
c        ----- calculate the auxiliary basis two-electron integrals
c
c calculate some goofy arrays for twoel
         call twoegoof(a(fna0),a(tmpsiz),a(lval),a(mycart),a(yn0),
     $        mmax,mindx)
         call setcc
         call twoel(z(c),z(exaux),z(rneed),a(need),z(xcont),z(s),
     $        a(xptprm),a(xnoprm),a(xnocon),a(xptcon),nat,npaux,
     $        maxcor-rneed+1,ntypes,nbtype,nnpx,ncontx,a(xstart),
     $        nbasx,a(nocart),a(nobf),a(maxmom),a(mintyp),
     $        a(nx),a(ny),a(nz),a(minmom),ops,mmax,a(fna0),
     $        a(tmpsiz),a(lval),a(mycart),a(yn0),mindx)
         if(logkey(ops,'print=int=auxtwo',.false.,' ')) then
            write(iout,120)
            call trtosq(z(t1),z(s),nbasx,nnpx)
            call wlmat(z(t1),nbasx,nbasx,xbflbl,xbflbl)
         end if
c 
c        ----- project the auxiliary two-electron integrals -----
         call trtosq(z(u),z(sinvx),nbasx,nnpx)
         call trtosq(z(t1),z(s),nbasx,nnpx)
         call ebc(z(t2),z(t1),z(u),nbasx,nbasx,nbasx)
         call ebtc(z(t1),z(u),z(t2),nbasx,nbasx,nbasx)
         write(iout,*)"triangularizing projected 2-e's..."
         call sqtotr(z(s),z(t1),nbasx,nnpx)
         write(iout,*)"...done"
         if(logkey(ops,'print=int=pauxtwo',.false.,' ')) then
            write(iout,130)
            call trtosq(z(u),z(s),nbasx,nnpx)
            call wlmat(z(u),nbasx,nbasx,xbflbl,xbflbl)
         end if
c
c        ----- send projected two-electron ints to disk -----
         call iosys('write real "auxiliary two-electron integrals"'
     $              //'to rwf',nnpx,z(s),0,' ')
      endif
c
c
      return
      end
