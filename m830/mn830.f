*deck @(#)mn830.f	5.1  11/6/94
      subroutine mn830(a,z,maxcor,tuniti,tunito,calc,
     $                 tfile,tfile1,tfile2)
c***begin prologue     mn830
c***date written       871027  (yymmdd)  
c***revision date      910708  (yymmdd)
c
c    8 july            rlm at lanl
c       adding separate input and output files.
c       fixing bug in the reordr if calc.eq.'ci'.
c***keywords           
c***author             saxe,paul (lanl) 
c***source             @(#)mn830.f	5.1   11/6/94
c***purpose            
c***description
c     
c    
c
c***references
c
c***routines called
c
c***end prologue       mn830
c
      implicit integer (a-z)
c
      character*(*) calc
      character*(*) tfile
      character*(*) tfile1
      character*(*) tfile2
      character*(*) tuniti,tunito
      character*16 inunit,ounit
      character*64 infile
      character*64 file1
      character*64 file2
      character*4096 ops
      integer a(maxcor)
      real*8 z(*)
c
      common /io/ inp,iout
c
      parameter (minbin=1024)
      parameter (minbuf=30000)
c
c
      file1=tfile1
      file2=tfile2
      inunit=tuniti
      ounit=tunito
      infile=tfile
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
      imngrp=1
      imxgrp=imngrp+ngroup
      jmngrp=imxgrp+ngroup
      jmxgrp=jmngrp+ngroup
      kadd=jmxgrp+ngroup
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
c     ----- fix up the symmetry array and add group mins and maxs -----
c
      call fixdrt(a(orbsym),norbs,a(ijgrp),nnp,a(imngrp),a(imxgrp),
     $     a(jmngrp),a(jmxgrp),ngroup)
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
      ints=h+n1int
      val=ints+nmax
      lenbin=minbin
      lab=wpadti(val+lenbin)
      bin=lab+lenbin
      rsort=iadtwp(bin+lenbin)
      isort=wpadti(rsort)
      lnsort=min(wptoin(n2int),maxcor-isort+1)
      need=isort+lnsort
c
      call getscm(need,a,maxcor,'sort',0)
c
c     ----- if a ci calculation, the orbitals have been reordered
c
c..rlm  is this ever correct?
c     if (calc.eq.'ci') then
c        do 2 i=1,norbs
c           a(orbtbf+i-1)=i
c2       continue
c     end if
c..rlm
c
      call tocan(z(h),n1int,z(ints),z(val),a(lab),a(bin),lenbin,
     #     z(rsort),a(isort),lnsort,
     #     a(isym),a(jsym),a(ijsym),nsym,nnpsym,
     #     a(ijklpt),nnqsym,a(nso),a(symoff),
     #     a(kadd),a(ladd),a(ijgrp),a(ijadd),
     #     norbs,numij,a(bftorb),nbf,nmax,ngroup,ops,
     #     a(ijpt),inunit,ounit,a(imngrp),a(imxgrp),a(jmngrp),
     $     a(jmxgrp),n2int,infile,file2,
     $     file1,a(orbsym),a(orbtbf))
c
c$$$      nnp=numij
c$$$      call iosys('read real "mo 1pdm" from ints',nnp,z,0,' ')
c$$$      write (iout,53)
c$$$ 53   format (' reordered mo one-electron integrals')
c$$$      call print(z,nnp,norbs,iout)
c$$$c
c$$$      call iosys('read real "mo one-electron integrals" from ints',
c$$$     $     nnp,z,0,' ')
c$$$      write (iout,54)
c$$$ 54   format (' original mo one-electron integrals')
c$$$      call print(z,nnp,norbs,iout)
c$$$c
c$$$      call iosys('read real "mo 2pdm" from ints',nnp**2,z,0,' ')
c$$$      write (iout,55)
c$$$ 55   format (' reordered mo two-electron integrals')
c$$$      call matout(z,nnp,nnp,nnp,nnp,iout)
c$$$c
c$$$      call iosys('read real "mo two-electron integrals" from ints',
c$$$     $     nnp**2,z,0,' ')
c$$$      write (iout,56)
c$$$ 56   format (' original mo two-electron integrals')
c$$$      call matout(z,nnp,nnp,nnp,nnp,iout)
c
c
      return
      end
