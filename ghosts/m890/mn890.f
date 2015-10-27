*deck @(#)mn890.f	1.1  11/30/90
      subroutine mn890(z,a,maxcor,tunit,calc,type)
c
      implicit integer (a-z)
c
      character*(*) type
      character*(*) calc
      character*(*) tunit
      character*16 unit
      character*4096 ops
      character*80 title
      integer a(maxcor)
      real*8 z(*)
c
      common /io/ inp,iout
c
      parameter (minbin=1024)
      parameter (minbuf=30000)
c
c
      call iosys('read character options from rwf',-1,0,0,ops)
c
      unit=tunit
c
c     ----- read dimensions etc from the drt -----
c
      call iosys('read integer "number of drt functions" from rwf',
     $     1,nbf,0,' ')
      call iosys('read integer "symmetries in ci" from rwf',
     $     1,nsym,0,' ')
      call iosys('read integer norbs from rwf',1,norbs,0,' ')
      call iosys('read integer symorb from rwf',1,symorb,0,' ')
      call iosys('read integer numij from rwf',1,numij,0,' ')
      call iosys('read integer ngroup from rwf',1,ngroup,0,' ')
      call iosys('read integer nmax from rwf',1,nmax,0,' ')
c
      nnpsym=nsym*(nsym+1)/2
      nnqsym=nnpsym*(nnpsym+1)/2
      nnp=nbf*(nbf+1)/2
c
c     ----- allocate core for the drt and other arays -----
c
      kadd=1
      ladd=kadd+symorb
      ijadd=ladd+symorb
      ijgrp=ijadd+numij
      bfsym=ijgrp+numij
      orbtbf=bfsym+nbf
      bftorb=orbtbf+norbs
      orbsym=bftorb+nbf
      nso=orbsym+norbs
      ijpt=nso+nsym
      ijklpt=ijpt+numij
      isym=ijklpt+nnqsym
      jsym=isym+nnpsym
      ijsym=jsym+nnpsym
      symoff=ijsym+nnpsym
      need=symoff+nsym
c
      call getscm(need,a,maxcor,'drt arrays',0)
c
      call iosys('read integer kadd from rwf',symorb,a(kadd),0,' ')
      call iosys('read integer ladd from rwf',symorb,a(ladd),0,' ')
      call iosys('read integer ijadd from rwf',numij,a(ijadd),0,' ')
      call iosys('read integer ijgrp from rwf',numij,a(ijgrp),0,' ')
      call iosys('read integer iout from rwf',nbf,a(bftorb),0,' ')
      call iosys('read integer orbsym from rwf',norbs,a(orbsym),0,' ')
      call iosys('read integer bfsym from rwf',nbf,a(bfsym),0,' ')
      call iosys('read integer orbtbf from rwf',norbs,a(orbtbf),0,' ')
c
c     ----- retrieve parameters we'll need -----
c
c      call iosys('read integer "number of symmetries" from rwf',1,junk,
c     #             0,' ')
c      if (junk.ne.nsym) then
c         call lnkerr('inconsistency in the number of symmetries')
c      end if
c      call iosys('read integer "number 1 ints" from rwf',1,n1int,0,' ')
c      call iosys('read integer "number 2 ints" from rwf',1,n2int,0,' ')
c
c
c      call iosys('read integer "number of so" from rwf',nsym,
c     #            a(nso),0,' ')
c      call iosys('read integer "symmetry pointer" from rwf',nnqsym,
c     #            a(ijklpt),0,' ')
c      call iosys('read integer "pair i symmetry" from rwf',nnpsym,
c     #            a(isym),0,' ')
c      call iosys('read integer "pair j symmetry" from rwf',nnpsym,
c     #            a(jsym),0,' ')
c      call iosys('read integer "pair symmetry" from rwf',nnpsym,
c     #            a(ijsym),0,' ')
c      call iosys('read integer "pair pointer" from rwf',numij,
c     #            a(ijpt),0,' ')
c
      a(nso)=norbs
      call sympt(nsym,nnpsym,nnqsym,a(isym),a(jsym),a(ijsym),a(ijpt),
     $     a(ijklpt),a(symoff),a(nso),n1int,n2int,nunij,norbs)
c
c     ----- work out some sizes we need for scratch arrays -----
c
      maxblk=0
      do 1 i=1,nsym
         maxblk=max(maxblk,a(nso+i-1))
    1 continue
c
c     ----- allocate core for sorting the integrals -----
c
      call getscm(0,a,maxcor,'?',0)
c
      h=wpadti(need)
      in=h+n1int
      lenbuf=max(minbuf,maxblk**2)
      lenbuf=min(lenbuf,maxblk**4)
      val=in+lenbuf
      lenbin=minbin
      lab=wpadti(val+lenbin)
      bin=lab+lenbin
      rsort=iadtwp(bin+lenbin)
      isort=wpadti(rsort)
      lnsort=min(wptoin(ngroup*nmax),maxcor-isort)
      need=isort+lnsort
c
      call getscm(need,a,maxcor,'sort',0)
c
c     ----- if a ci calculation, the orbitals have been reordered
c
      write(iout,*) calc
      write(iout,*) (a(bftorb+i-1),i=1,nbf)
      if (calc.eq.'ci') then
         do 2 i=1,nbf
            a(bftorb+i-1)=i
 2       continue
      end if
c
c
c     open an output file for the integrals.
      intout=50
      open(unit=intout,file='intout',access='sequential',
     $     form='formatted',err=99,status='new')
      title='cu2o7..madelung potential..ecp background..dz-ecp2/3s2p'
      write(intout,1000) title
      write(intout,1001) norbs
 1000 format(a80)
 1001 format(10i5)
c
c
      call toguga(z(h),n1int,z(in),lenbuf,z(val),a(lab),a(bin),lenbin,
     #            z(rsort),a(isort),lnsort,
     #            a(isym),a(jsym),a(ijsym),nsym,nnpsym,
     #            a(ijklpt),nnqsym,a(nso),a(symoff),
     #            a(kadd),a(ladd),a(ijgrp),a(ijadd),
     #            norbs,numij,a(bftorb),nbf,nmax,ngroup,ops,intout,
     #            a(ijpt),unit,type)
c
      close(unit=intout)
      return
c
   99 call lnkerr('cannot open the output integral file.')
c
      return
      end
