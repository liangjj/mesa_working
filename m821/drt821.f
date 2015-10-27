*deck @(#)drt821.f	5.1  11/6/94
      subroutine drt821(    kadd,ladd,ijadd,ijgrp
     #,                        ningrp,orbsym,nblkmn,nblkmx
     #,                 ijxx,klxx,nklxx,ijww,klww,nklww,dsk)
c
      implicit integer (a-z)
c
      character*(*) dsk
      real*8 ver4x,drtver
c
      common /lblm/ lblint(26),lbldrt(26),ver4x,drtver
      common /io/     inp,iout
      common /dimn/ nbf,nsym,norbs,nrows,nrows4,nwks,nwks2,nlevs
     *,             nrowoc,nrow4o,nwksoc,nlevoc,levfrm
     *,             nwksmx,nlwkmx,nuwkmx,maxb,nroots,lvfrm1,nrefs
      common /intm/ nmax,ngroup,nblkoc,numij,symorb,intsrt
      common /x4x821/ nijvir
c
      dimension ijadd(numij),kadd(symorb),ladd(symorb)
      dimension ijgrp(numij)
      dimension orbsym(norbs)
      dimension ningrp(ngroup),nblkmn(ngroup),nblkmx(ngroup)
      dimension ijxx(numij ),klxx(nijvir),nklxx(nsym)
      dimension ijww(numij ),klww(nijvir),nklww(nsym)
c
      call iosys('read integer kadd from '//dsk,-1,kadd,0,' ')
      call iosys('read integer ladd from '//dsk,-1,ladd,0,' ')
      call iosys('read integer ijadd from '//dsk,-1,ijadd,0,' ')
      call iosys('read integer ijgrp from '//dsk,-1,ijgrp,0,' ')
      call iosys('read integer ningrp from '//dsk,-1,ningrp,0,' ')
      call iosys('read integer orbsym from '//dsk,-1,orbsym,0,' ')
      call iosys('read integer ijxx from '//dsk,-1,ijxx,0,' ')
      call iosys('read integer klxx from '//dsk,-1,klxx,0,' ')
      call iosys('read integer nklxx from '//dsk,nsym,nklxx,0,' ')
      call iosys('read integer ijww from '//dsk,-1,ijww,0,' ')
      call iosys('read integer klww from '//dsk,-1,klww,0,' ')
      call iosys('read integer nklww from '//dsk,nsym,nklww,0,' ')
c
c     ----- add offset for group to ijww and ijww -----
c
      do 30 i=1,numij
         junk=(ijgrp(i)-1)*nmax
         ijww(i)=ijww(i)+junk
         ijxx(i)=ijxx(i)+junk
   30 continue
c
c
c
c
c     write (iout,*) ' kadd '
c     write (iout,710) (kadd(i),i=1,symorb)
c     write (iout,*) ' ladd '
c     write (iout,710) (ladd(i),i=1,symorb)
c     write (iout,*) ' ijadd'
c     write (iout,710) (ijadd(i),i=1,numij)
c     write (iout,*) ' ijgrp'
c     write (iout,710) (ijgrp(i),i=1,numij)
c     write (iout,*) ' ningrp'
c     write (iout,710) (ningrp(i),i=1,ngroup)
c     write (iout,*) ' orbsym'
c     write (iout,710) (orbsym(i),i=1,norbs)
c     write (iout,*) ' ijxx '
c     write (iout,710) (ijxx(i),i=1,numij)
c     write (iout,*) ' klxx '
c     write (iout,710) (klxx(i),i=1,nijvir)
c     write (iout,*) ' nklxx'
c     write (iout,710) (nklxx(i),i=1,nsym)
c     write (iout,*) ' ijww '
c     write (iout,710) (ijww(i),i=1,numij)
c     write (iout,*) ' klww '
c     write (iout,710) (klww(i),i=1,nijvir)
c     write (iout,*) ' nklww'
c     write (iout,710) (nklww(i),i=1,nsym)
c
c
c
  710 format (20i5)
    9 continue
      iblk=1
      return
      end
