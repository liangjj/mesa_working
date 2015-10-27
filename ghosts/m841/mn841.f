*deck @(#)mn841.f	1.1  11/30/90
      subroutine mn841(a,z,maxcor)
c
      implicit integer (a-z)
c
      character*4096 ops
      integer a(maxcor)
      logical mcscf
      logical ci
      logical logkey
      real*8 z(*)
c
      common /io/ inp,iout
c
      parameter (minbuf=30000)
c
c
      mcscf=logkey(ops,'mcscf',.false.,' ')
      ci=logkey(ops,'ci',.false.,' ')
c
c     ----- read dimensions etc -----
c
      call iosys('read integer "number of basis functions" from rwf',
     $     1,nbf,0,' ')
      call iosys('read integer "number of grouped 2pdm elements" '//
     $     'from rwf',1,n2pdm,0,' ')
      call iosys('read integer "number of momentum groups" from rwf',
     $     1,ngrp,0,' ')
c
      nnp=nbf*(nbf+1)/2
      nnpgrp=ngrp*(ngrp+1)/2
c
c     ----- allocate core for the pointer arrays -----
c
      ntriang=min(nnp,minbuf/nnp+1)
c
      bftgrp=1
      bftcmp=bftgrp+nbf
      grpsiz=bftcmp+nbf
      gptij=grpsiz+ngrp
      gptkl=gptij+nnpgrp
      gklsiz=gptkl+nnpgrp
      bin=gklsiz+nnpgrp
      lab=bin+nnp
      val=iadtwp(lab+nnp)
      dm=val+nnp
      rsort=dm+nnp*ntriang
      isort=wpadti(rsort)
c
      call getscm(0,a,maxcor,'?',0)
c
      lnsort=min(wptoin(n2pdm),maxcor-isort)
c
      need=isort+lnsort
c
      call getscm(need,a,maxcor,'pointer arrays',0)
c
      call iosys('read integer "basis functions to groups" from rwf',
     $     nbf,a(bftgrp),0,' ')
      call iosys('read integer "basis functions to components" '//
     $     'from rwf',nbf,a(bftcmp),0,' ')
      call iosys('read integer "momentum group size" from rwf',
     $     ngrp,a(grpsiz),0,' ')
      call iosys('read integer "group ij pointers" from rwf',
     $     nnpgrp,a(gptij),0,' ')
      call iosys('read integer "group kl pointers" from rwf',
     $     nnpgrp,a(gptkl),0,' ')
      call iosys('read integer "group pair blocksize" from rwf',
     $     nnpgrp,a(gklsiz),0,' ')
c
c
      call toder(z(dm),z(val),a(lab),a(bin),z(rsort),a(isort),lnsort,
     $     a(bftgrp),a(bftcmp),a(grpsiz),a(gptij),a(gptkl),a(gklsiz),
     $     nbf,nnp,ngrp,nnpgrp,ntriang,n2pdm,ci,mcscf)
c
      write(iout,99011)
99011 format(/,5x,'  m841: group ordered 2e-ints completed')
c
      return
      end
