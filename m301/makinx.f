*deck  @(#)makinx.f	1.3 7/30/91
      subroutine makinx(pkindx,nbf,newnbf,newnnp,bflabl,nulabl,
     $                  prnt)
      implicit integer (a-z)
      real*8 z
      integer a
      pointer(p,z(1)), (p,a(1))
      dimension pkindx(nbf)
      character*16 bflabl(*),nulabl(*)
      logical prnt
      common /io/ inp,iout
c
      call iosys('read integer "number of atoms" from rwf',1,nat,0,' ')
      call iosys('length of exponents on rwf',nprim,0,0,' ')
      call iosys('length of "contraction coefficients" on rwf',
     $     ncont,0,0,' ')
      call iosys('length of "number of pure functions" on rwf',
     $     ntypes,0,0,' ')
      call iosys('read integer "number basis types" from rwf',
     $     1,nbtype,0,' ')
      call iosys('length of "power of x" on rwf',lenxyz,0,0,' ')
c
      nnp=(nbf+1)*nbf/2
c
c     divide core for basis set information.
c
      zan=1
      c=zan+nat
      cont=c+3*nat
      ex=cont+ncont
      ptprim=wpadti(ex+nprim)
      noprim=ptprim+ntypes*nat
      nocont=noprim+ntypes*nat
      ptcont=nocont+ntypes*nat
      start=ptcont+ntypes*nat
      nocart=start+ntypes*nat
      bstart=nocart+ntypes*nat
      blen=bstart+ntypes*nat
      nstart=blen+ntypes*nat
      lenbf=nstart+nat
      nobf=lenbf+nat
      minmom=nobf+ntypes
      maxmom=minmom+ntypes
      mintyp=maxmom+ntypes
      nx=mintyp+ntypes
      ny=nx+lenxyz
      nz=ny+lenxyz
      need=nz+lenxyz
c
c     make sure we have this.
c      call getscm(need,a,maxcor,'makinx',0)
c      call manmem(need,p,ngot,'makinx',0)
      call getmem(need,p,ngot,'makinx',0)
c
c     determine the amount of scratch space required.
c     retrieve information about the most demanding shell block.
c
      call iosys('read integer maxprm from rwf',1,maxprm,0,' ')
      call iosys('read integer maxcont from rwf',1,mxcont,0,' ')
      call iosys('read integer maxl from rwf',1,maxl,0,' ')
      call iosys('read integer maxblk from rwf',1,maxblk,0,' ')
c
c     read in basis set information from read-write file.
c
      call iosys('read real exponents from rwf',-1,z(ex),0,' ')
      call iosys('read real "contraction coefficients" from rwf',
     $     -1,z(cont),0,' ')
      call iosys('read real "nuclear charges" from rwf',-1,z(zan),0,' ')
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
      call iosys('read integer "pointer to cartesians" from rwf',
     $     -1,a(mintyp),0,' ')
      call iosys('read integer "pointer to first function" from rwf',
     $     -1,a(start),0,' ')
      call iosys('read integer "power of x" from rwf',-1,a(nx),0,' ')
      call iosys('read integer "power of y" from rwf',-1,a(ny),0,' ')
      call iosys('read integer "power of z" from rwf',-1,a(nz),0,' ')
c
c
c
      call iosys('read character "basis function labels" from rwf',
     #            -1,0,0,bflabl)
c
      call basout(z(c),z(ex),z(cont),a(ptprim),
     #            a(noprim),a(nocont),a(ptcont),nat,nprim,
     #            ntypes,nbtype,nnp,ncont,a(start),
     #            z(zan),a(nocart),a(nobf),a(maxmom),a(mintyp),
     #            a(nx),a(ny),a(nz),a(minmom),a(nstart),a(bstart),
     $            a(blen),bflabl,prnt)
c
      call setind(a,pkindx,nbf,newnbf,bflabl,nulabl,
     #           a(nstart),a(bstart),a(blen),a(nocart),a(nocont),nat,
     #           ntypes,prnt,a(lenbf))
c
      newnnp=newnbf*(newnbf+1)/2
c
c      call manmem(-ngot,p,idum,'makinx',idum)
      call getmem(-ngot,p,idum,'makinx',idum)
      return
      end
