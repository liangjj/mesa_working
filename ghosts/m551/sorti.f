*deck %W%  %G%
      subroutine sorti(z,a,nwint,g1,g2,intape)
c
c***begin prologue
c***date written       871022   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c
c***keywords
c***author             lengsfield, byron (brl)
c***source             %W%   %G%
c
c***purpose
c
c***description
c
c***references
c
c***routines called    (none)
c
c***end prologue
c
c
c
      implicit integer (a-z)
c
      real*8 z,g1,g2
      integer a(1)
c...
c...  equivalence (z,a)
c...  integer a(nwint)
c...
c
      dimension z(2),g1(2),g2(2)
c
      common /io/ inp,iout
c
c
c     ----- read dimensions etc from drt tape -----
c
      call iosys('read integer "number of drt functions" from rwf',
     $     1,nbfs,0,' ')
      call iosys('read integer "symmetries in ci" from rwf',
     $     1,nsym,0,' ')
      call iosys('read integer norbs from rwf',1,norbs,0,' ')
      call iosys('read integer nrows from rwf',1,nrows,0,' ')
      nrows4=nrows*4
      call iosys('read integer nlevs from rwf',1,nlevs,0,' ')
      call iosys('read integer nrefs from rwf',1,nrefs,0,' ')
      call iosys('read integer orbfrm from rwf',1,levfrm,0,' ')
      call iosys('read integer symorb from rwf',1,symorb,0,' ')
      call iosys('read integer numij from rwf',1,numij,0,' ')
      call iosys('read integer ngroup from rwf',1,ngroup,0,' ')
      call iosys('read integer nmax from rwf',1,nmax,0,' ')
      call iosys('read integer nijvir from rwf',1,nijvir,0,' ')
c
cps
c
c     ----- divide core and get the drt arrays -----
c
      kadd=27
      ladd=kadd+symorb
      ijadd=ladd+symorb
      ijgrp=ijadd+numij
      bftorb=ijgrp+numij
      orbsym=bftorb+nbfs
      ijww=orbsym+norbs
      klww=ijww+numij
      ijxx=klww+nijvir
      klxx=ijxx+numij
      need=klxx+nijvir
c
c     ----- dynamic core allocation by pws -----
c
      call getscm(need,a,maxcor,'sorti',0)
c
c
      call mctdrt(a(kadd),a(ladd),a(ijadd),a(ijgrp),
     $     a(bftorb),a(orbsym),end,ngroup,nrefs,symorb,numij,nbfs,
     $     norbs,itape8,nijvir,a(ijww),a(klww),a(ijxx),a(klxx),
     $     levfrm,nlevs,nrows,nsym)
c
c        write(iout,*)' nbf nbfs ',nbf,nbfs
c
      nbf=nbfs
c
      hmo=iadtwp(need)
      c=hmo+numij
c
      nnp=(nbf+1)*nbf/2
c
      nnb=norbs*(norbs+1)/2
      t1=c+nbf*nbf
      t2=t1+nbf**2
      values=t2+nnb*nbf*nbf
c
      ntriang=nnb
c
c
c       temporarily get integrals off of tape 99
c      later  they will be passed in
c
c
c
c
c      val=values+ntriang*nnp
c
      val=values
      lab=wpadti(val+ntriang*nnb+2)
      bin=iadtwp(lab+ntriang*nnb+2)
c
      maxipt=ntriang*nnb+2
c==================================================================c
cbhl
      call getscm(0,a,maxcor,'ssorti',0)
      nwwp=iadtwp(maxcor)-512
cbhl
      lenbin=ntriang*nnp
      asort=bin+lenbin
      lnsort=max(nwwp-asort,2000)
      lnsort=max(lnsort,10000)
c
c     ----- ps trying to limit core size -----
c
cps      lnsort=min(lnsort,nnp**2)
      need=asort+lnsort
c
c     ----- dynamic core allocation by pws -----
c
      call getscm(wpadti(need),a,maxcor,'sorti',0)
c
cmp
cmp
c
c.io     itap50=intape
cps      call sfile(itap50,'tape50')
cps      call ncdlbl(itap50,a,a,ngroup,nmax,nsym,0.0d+00,ecore,enuc)
cps      call getwad(itap50,st50)
c******************
c.io      itape2=itap50
c.io      st2=st50
c******************
c.io      itap91=91
cps      call sfile(itap91,'tape91')
cps      call srew(itap91)
cps      call getwad(itap91,st91)
c
c      write(iout,*)'  calling srtdrt'
c c c  write(iout,*)' values c t1 t2 val lab asort bin'
c c c  write(iout,*) values,c,t1,t2
c c c  write(iout,*) val,lab,asort,bin
c
      call srtdrt(g2,nnp,ntriang,z(c),nbfs,norbs,z(t1),z(t2),
     $     numij,ngroup,nmax,nsym,a(ijgrp),a(ijadd),a(kadd),a(ladd),
     $     a(orbsym),g1,levfrm,a(ijww),a(klww),a(ijxx),
     $     a(klxx),nijvir, z(val),a(lab),maxipt,
     $     z(asort),lnsort,itape2,st2,itap91,st91,lenbin,z(bin))
c
c     write(iout,9091)
 9091 format(/,'  srtdrt completed ')
c
c
ccccccmaxi=maxcor*intowp(1)
cccccccallczero(z,need)
ccccccwrite(iout,*)c'cccisortccalledc'
cccccccallccisort(z,a,maxi,itap50)
ccccccwrite(iout,*)c'cccimaincccalledc'
cccccciflag=0
cccccccallccimain(z,a,maxcor,iflag,z,z,itap50)
ccccccwrite(iout,*)c'ccendcccimainc'
c
c     ----- timing -----
c
c      call endtim
c
c
      return
      end
