*deck @(#)mn821.f	5.2  2/5/95
      subroutine mn821(file,unit)
c***begin prologue     mn821.f
c***date written       840823     yymmdd  
c***revision date      7/27/92      
c
c  25 february  1991    rlm at lanl
c     removing limitations on maxcor imposed below.
c
c  18 march 1988   bhl at llnl
c        maxpos limited to 2,000,000 for cos io
c
c  29 april 1985        pws at lanl
c        mesax version implemented on cray x-mp
c
c
c***keywords          sort,guga,integrals 
c***author            saxe, paul(lanl)
c***source             @(#)mn821.f	5.2   2/5/95
c***purpose: to create the 3 & 4-external element supermatrices
c            for the guga matrix ci program. this program reads
c            the transformed integrals from tape2, forms the
c            supermatrices and places them in the reserved spaces
c            on tape2.
c***description
c     
c    
c
c***references
c
c***routines called
c
c***end prologue       mn821.f
c
      implicit integer (a-z)
c
      character*(*) file, unit
      character*4096 ops
      character*8 prtflg, dsk
      character*16 chrkey, key
      character*128 nmcnfg
      logical prnt
      real*8 z
      real*8 rep,fzcore,eguess,eci,cnverg,sqcdif
      real*8 drtver,ver4x
      integer ia, iz
      integer ngot(2)
      pointer(p,ia(1))
      pointer(pz,z(1)), (pz,iz(1))
c
      common /io/      inp,iout
      common /dimn/   nbf,nsym,norbs,nrows,nrows4,nwks,nwks2,nlevs
     *,               nrowoc,nrow4o,nwksoc,nlevoc,levfrm
     *,               nwksmx,nlwkmx,nuwkmx,maxb,nroots,lvfrm1,nrefs
      common /intm/   nmax,ngroup,nblkoc,numij,symorb,intsrt
      common /diag/   rep,fzcore,eguess,eci,refwlk,mxiter,cnverg,icnvg
     *,               iter,sqcdif,nroot
      common /lblm/   lblint(26),lbldrt(26),ver4x,drtver
      common /x4x821/   nijvir
      common /sort/   iblock,lnbuf,maxsrt
c
      equivalence (drtver,iver)
c
      data prnt/.true./
      save prnt
c
 1000 format(1x,'shape driven guga sort:')
c
c     ----- recover the options string -----
c
      write(iout,*) ' m821: Guga Sort'
      call iosys('read character options from rwf',-1,0,0,ops)
      key=chrkey(ops,'int=drt=key','drt',' ')
      call pakstr(key,lenkey)
      call iosys('read character "drt file name '//key(1:lenkey)
     1            //'" from rwf',0,0,0,dsk)
      call iosys('read character "hamiltonian manipulation '//
     $           'filename" from rwf',0,0,0,nmcnfg)
      call iosys('open hconfig as old',0,0,0,nmcnfg)
      call getmem(0,p,ngot,'first',0)
      call iosys('read integer mncore from rwf',1,bottom,0,' ')
      call iosys('read integer mxcore from rwf',1,top,0,' ')
      write(iout,*) 'minimum core available = ',bottom
      write(iout,*) 'maximum core available = ',top
      if(chrkey(ops,'core-use=m821','minimum',' ').eq.'minimum') then
         top=bottom
      endif
c      call getmem(top,p,maxcor,'m821',1)
c
c     ----- start the timing routines -----
c
c
c     has printing been turned off externally?
c
      maxi=top
      maxr=iadtwp(top)
      call iosys('read character "print flag" from rwf',-1,0,0,prtflg)
      if(prtflg.eq.'minimum') prnt=.false.
c
c
      ver4x=5.290485
c
      if(prnt) then
         write(iout,1000)
      endif
c
c     ----- read dimensions etc from read-write file -----
c
      call iosys('read integer "number of drt functions" '//
     $           'from '//dsk,1,nbf,0,' ')
      call iosys('read integer "symmetries in ci" from '//dsk,
     $            1,nsym,0,' ')
      call iosys('read integer norbs from '//dsk,1,norbs,0,' ')
      call iosys('read integer nrows from '//dsk,1,nrows,0,' ')
      call iosys('read integer nlevs from '//dsk,1,nlevs,0,' ')
      call iosys('read integer nrefs from '//dsk,1,nrefs,0,' ')
      call iosys('read integer orbfrm from '//dsk,1,lvfrm1,0,' ')
      call iosys('read integer symorb from '//dsk,1,symorb,0,' ')
      call iosys('read integer numij from '//dsk,1,numij,0,' ')
      call iosys('read integer ngroup from '//dsk,1,ngroup,0,' ')
      call iosys('read integer nmax from '//dsk,1,nmax,0,' ')
      call iosys('read integer nijvir from '//dsk,1,nijvir,0,' ')
      call iosys('read integer nwks from '//dsk,1,nwks,0,' ')
c
      levfrm=lvfrm1+1
      if (prnt) then
c      write(iout,31)lbldrt
   31 format(5x,'label from drt...',26a3)
      write(iout,32)nwks
   32 format(5x,'number of configurations:     ',i8)
      write(iout,33)nbf
   33 format(5x,'number of basis functions:    ',i8)
      write(iout,34) norbs
   34 format(5x,'number of allowed orbitals:   ',i8)
      write(iout,35)lvfrm1
   35 format(5x,'fermi level:                  ',i8)
      end if
c
c     ----- calculate pointers into blank common -----
c
      ijxx=1
      klxx=ijxx+numij
      nklxx=klxx+nijvir
      ijww=nklxx+nsym
      klww=ijww+numij
      nklww=klww+nijvir
      ijadd=nklww+norbs
      ijgrp=ijadd+numij
      kadd=ijgrp+numij
      ladd=kadd+symorb
      ningrp=ladd+symorb
      nblkmn=ningrp+ngroup
      nblkmx=nblkmn+ngroup
      orbsym=nblkmx+ngroup
      leni=orbsym+norbs
      need = leni + 10000
      call getmem(need,p,ngot(1),'m821 drt',1)
      maxi = maxi - ngot(1)
      maxr=iadtwp(maxi)
c
C      if (leni.gt.maxi) then
c        add a little to be on the safe side.
c      call getscm(leni,z,maxi,'m321 getdrt',0)
c         maxr=iadtwp(maxi)
c      end if
c
      call drt821(        ia(kadd),ia(ladd),ia(ijadd),ia(ijgrp)
     #,                                                  ia(ningrp)
     #,        ia(orbsym),ia(nblkmn),ia(nblkmx)
     #,        ia(ijxx),ia(klxx),ia(nklxx),ia(ijww),ia(klww),
     #         ia(nklww),dsk)
c
c      ints=need
c
      ints=1
      h3=ints+nmax
c      asort=h3+nijvir
      val=h3+nijvir
      lenbin=min(ngroup*nmax,2000000)
      lab=wpadti(val+lenbin)
      bin=lab+lenbin
      offset=bin+lenbin
      minsym=offset+norbs
      maxsym=minsym+nsym
      isort=maxsym+nsym
      asort=iadtwp(isort)
      need = asort + ngroup*nmax + 1000
      if(need.gt.maxr) then
         need=wpadti(maxcor)
         call getmem(need,pz,ngot(2),'m821 srt',1)
      else
         need=wpadti(need)
         call getmem(need,pz,ngot(2),'m821 srt',1)
      endif
      lnsort=ngot(2)
      write(iout,1) ngroup*nmax, lenbin, lnsort
 1    format(/,1x,'amount of data to sort = ', i10,
     1       /,1x,'length of sorting bin  = ',i10,
     2       /,1x,'free core for sort     = ',i10)
c
c      if (asort+ngroup*nmax.gt.maxr) then
c         call getscm(0,z,maxpos,'mn821 how much?')
c.bhl.cos.io
c         maxpos=iadtwp(maxpos)-1
c        maxpos=min(maxpos,2000000)
c.bhl.cos.io
c        add just a little to be on the safe side
c         need=min(maxpos,asort+ngroup*nmax+100)
c         call getscm(wpadti(need),z,maxr,'mn821 sort',0)
c         maxr=iadtwp(maxr)
c      end if
c
      
c
c     ----- initialize the sorting routine -----
c
      call sorter('start',z(asort),iz(isort),lnsort,ngroup*nmax,
     #             1024,-1,0,0,file,unit,prnt)
c
      call xxint(z(ints),z(h3),ia(orbsym),ia(ijgrp),ia(ijadd)
     #,          ia(kadd),ia(ladd),ia(ijww),ia(klww),ia(ijxx),
     #           ia(klxx),z(val),iz(lab),iz(bin),iz(offset),
     #           iz(minsym),iz(maxsym),lenbin,z(asort),prnt)
c
      call sorter('end',z(asort),iz(isort),0,0,0,0,0,0,0,0,prnt)

c
c     ----- end timing -----
c
c
c     ----- and exit with dignity -----
c
      call getmem(-ngot(1),p,idum,'m821 drt',idum)
      call getmem(-ngot(2),pz,idum,'m821 srt',idum)
      return
  920 write(iout,925)norbs,norb2
  925 format(' number of orbitals does not match ',2i10)
      call lnkerr(' ')
  930 write(iout,935)leni
  935 format(' insufficient space in integer blank common need',i10,
     *' words.')
      call lnkerr(' ')
  940 write(iout,945)lenr
  945 format(' insufficient space in blank common need',i10,' words.')
      end
