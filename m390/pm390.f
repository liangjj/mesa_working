*deck @(#)pm390.f	1.1  8/1/91
      subroutine pm390(z,a)
c
      implicit integer (a-z)
c
      parameter (maxsym=8,maxao=500)
c
      real*8 z(1)
      real*8 nucrep
      integer blabel(20),nd(maxsym),ityp(maxsym),nso(maxsym)
      integer mtype(2,maxsym),ms(142),mnl(142)
      integer bfsym(maxao),bfnum(maxao)
      character*4096 ops
      character*128 namint
      character *8 chrkey,cpass, ints14, gvbvec
      dimension ititle(18),ix(10)
      integer orbsym(maxao)
      character *800 card
c
      common /io/ inp,iout
      common /lbf/ kbuff(20000b) , jbuff(2000b) , ibuff(10000b)
c     common a(1)
c
c     ----- work out the symmetry pointer arrays -----
c
      call posinp('$symtrn',cpass)
      call cardin(card)
      ints14=chrkey(card,'ints14','ints14',' ')
      gvbvec=chrkey(card,'gvbvec','gvbvec',' ')
      nsym=intkey(card,'no-symmetries',1,' ')
c
      write (iout,11)
   11 format(/,10x,'m810:symmetry sort and scf blocking-c2v symmetry')
      write(iout,2) ints14,gvbvec
    2 format(/,5x,'input integrals file:',1x,a8,1x,'transformation matri
     1x:',1x,a8)
c     ----- read header information from integral tape -----
c*****   integral tape is from the canonical sorting of ijkl
c*****
c
      call qassign(14,ints14,kbuff,20000b)
      rewind 14
      read(14) ititle
      read(14) nao
      read(14)
      read(14)
      read(14) nucrep
      read(14) (ix(i),i=1,10),k
      nvec=2*k
      do 60 i=1,nvec
         read(14)
 60   continue
      write(iout,25) nao,nsym,nucrep
   25 format(/,5x,'number atomic orbs',1x,i4,1x,'number symmetries',1x,
     1  i2,1x,/,5x,'nuclear repulsion',1x,e13.6)
c
c     ----- how many basis functions are there? -----
c
      nbf=nao
      do 65 i=1,nbf
         orbsym(i)=1
         bfsym(i)=1
         bfnum(i)=i
 65   continue
      call intarr(card,'orbital-sym',orbsym,nao,' ')
      call intarr(card,'basis-fn-sym',bfsym,nao,' ')
c
c
c 
c ***** bfsym(i)=the symmetry number of the ith ao
c ***** bfnum(i)= the number of the ith ao in it's symmetry group
      write(iout,3) (i,bfsym(i),i=1,nao)
    3 format(/,5x,'atomic symmetries',(/,5x,10('(',i2,','i2,')':)))
      write(iout,4) (i,orbsym(i),i=1,nao)
    4 format(/,5x,'molecular symmetries',(/,5x,10('(',i2,','i2,')':)))
c 
      call izero(nso,nsym)
      do 66 i=1,nbf
         j=bfsym(i)
         nso(j)=nso(j)+1
         bfnum(i)=nso(j)
 66   continue
c 
      write(iout,29) (bfnum(i),i=1,nao)
   29 format(/,5x,'atomic orbital number',( /,10(1x,i4,1x)))
c 
      nnp=nbf*(nbf+1)/2
c
c     ----- divvy up core for symmetry pointer arrays -----
c
      nnpsym=nsym*(nsym+1)/2
      nnqsym=nnpsym*(nnpsym+1)/2
c
      isym=1
      jsym=isym+nnpsym
      ijsym=jsym+nnpsym
      ijpt=ijsym+nnpsym
      ijklpt=ijpt+nnp
      symoff=ijklpt+nnqsym
      need=symoff+nsym
c
      call getscm(need,a,maxcor,' ',0)
c
      call sympt(nsym,nnpsym,nnqsym,a(isym),a(jsym),a(ijsym),a(ijpt),
     #     a(ijklpt),a(symoff),nso,n1int,n2int,nnp,nbf)
c
c     ----- create an iosys integral file -----
c
      call iosys('read character "integral filename" from rwf',
     #     -1,0,0,namint)
      call iosys('open ints as new',0,0,0,namint)
c
c     ----- store these pointer and other useful arrays -----
c
      call iosys('write real "nuclear repulsion energy" to rwf',
     $     1,nucrep,0,' ')
      call iosys('write integer "number of basis functions" to rwf',1,
     $     nbf,0,' ')
      call iosys('write integer "number of so" to rwf',nsym,nso,0,' ')
      call iosys('write integer "number of symmetries" to rwf',1,nsym,
     #     0,' ')
      call iosys('write integer "number 1 ints" to rwf',1,n1int,0,' ')
      call iosys('write integer "number 2 ints" to rwf',1,n2int,0,' ')
      call iosys('write integer "pair i symmetry" to rwf',nnpsym,
     #     a(isym),0,' ')
      call iosys('write integer "pair j symmetry" to rwf',nnpsym,
     #     a(jsym),0,' ')
      call iosys('write integer "pair symmetry" to rwf',nnpsym,
     #     a(ijsym),0,' ')
      call iosys('write integer "pair pointer" to rwf',nnp,a(ijpt),
     #     0,' ')
      call iosys('write integer "symmetry pointer" to rwf',nnqsym,
     #     a(ijklpt),0,' ')
      call iosys('write integer "symmetry offsets" to rwf',nsym,
     #     a(symoff),0,' ')
c
c ***** change lenbuf to 500 since integrals are read in , in 500 blocks
      lenbuf=500
c *****
c
c     ----- core for one-electron part -----
c
      t1=iadtwp(need)
      t2=t1+n1int
      rbuf=t2+n1int
      ibuf=wpadti(rbuf + lenbuf)
      need=ibuf+lenbuf
c
      call getscm(need,a,maxcor,' ',0)
c
c     ----- overlap integrals -----
c
      call oneint(z(t1),n1int,a(ibuf),z(rbuf),lenbuf,nnp,a(ijpt),
     #     a(symoff),nsym,itape,bfsym,bfnum)
      call iosys('write real "so overlap integrals" to rwf',
     $     nnp,z(t1),0,' ')
c
c      write (iout,200) 
c 200  format (5x,'the overlap integrals:')
c      call print(z(t1),nnp,nbf,iout)
c
c     ----- one-electron integrals -----
c
      call oneint(z(t1),n1int,a(ibuf),z(rbuf),lenbuf,nnp,a(ijpt),
     #     a(symoff),nsym,itape,bfsym,bfnum)
c
      call iosys('write real "so one-electron integrals" to rwf',
     $     n1int,z(t1),0,' ')
c
      call oneint(z(t2),n1int,a(ibuf),z(rbuf),lenbuf,nnp,a(ijpt),
     #     a(symoff),nsym,itape,bfsym,bfnum)
c
      call iosys('write real "sokinetic integrals" to rwf',
     $     n1int,z(t2),0,' ')
c
c      write (iout,202)
c 202  format(/,5x,'the kinetic integrals:')
c      call print(z(t2),nnp,nbf,iout)
c
c
      do 1 i=1,n1int
         z(t1+i-1)=z(t1+i-1)-z(t2+i-1)
    1 continue
c
      call iosys('write real "so potential integrals" to rwf',
     $     n1int,z(t1),0,' ')
c
c      write (iout,201)
c 201  format(/,5x,'the potential integrals:')
c      call print(z(t1),nnp,nbf,iout)
c
c     ----- core for two-electron part -----
c
      lenbin=1024
c
      rbuf=iadtwp(need)
      ibuf=wpadti(rbuf+lenbuf)
      values=iadtwp(ibuf+lenbuf)
      ptr=wpadti(values+lenbin)
      bins=ptr+lenbin
      sort=iadtwp(bins+lenbin)
      asort=wpadti(sort)
c
      call getscm(0,a,maxcor,' ',0)
c
      lnsort=min(maxcor-asort+1,n2int)
      need=asort+lnsort
c
      call getscm(need,a,maxcor,' ',0)
c
      call twoint(a(ibuf),z(rbuf),lenbuf,z(values),a(ptr),lenbin,
     #     z(sort),a(asort),lnsort,n2int,a(ijklpt),nnqsym,nso,
     #     nsym,itape,a(bins),intfil,bfsym,bfnum)
c
      call count(a(ibuf),z(rbuf),lenbuf,z(values),a(ptr),lenbin,
     #     z(sort),a(asort),lnsort,n2int,a(ijklpt),nnqsym,nso,
     #     nsym,itape,a(bins),intfil)
c
c     ----- symmetry block the vector -----
c
      nblock=0
      do 100 i=1,nsym
         nblock=nblock+nso(i)**2
 100  continue
c
      cao=1
      cso=cao+nbf**2
      eigval=cso+nblock
      need=wpadti(eigval+nbf)
c
      call getscm(need,a,maxcor,' ',0)
c
      call cblock(gvbvec,z(cao),z(cso),bfsym,bfnum,orbsym,nbf,nblock,nso
     $     ,nsym,ibuff,z(eigval))
c
      call chainx(0)
c
c
      return
      end
