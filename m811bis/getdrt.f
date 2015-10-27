*deck @(#)getdrt.f	1.1  11/30/90
      subroutine getdrt(kadd,ladd,ijadd,ijgrp,iout,orbsym
     #,                 end,ngroup,nrefs,symorb,numij,nbf,norbs,
     #                  levfrm,nlevs,nrows,nsym,dsk)
c
c  read in drt information
c
      implicit integer (a-z)
c
      integer kadd(symorb),ladd(symorb),ijadd(numij),ijgrp(numij)
      integer iout(nbf),orbsym(norbs)
      character*(*) dsk
c
      call iosys('read integer kadd from '//dsk,-1,kadd,0,' ')
      call iosys('read integer ladd from '//dsk,-1,ladd,0,' ')
      call iosys('read integer ijadd from '//dsk,-1,ijadd,0,' ')
      call iosys('read integer ijgrp from '//dsk,-1,ijgrp,0,' ')
      call iosys('read integer iout from '//dsk,-1,iout,0,' ')
      call iosys('read integer orbsym from '//dsk,-1,orbsym,0,' ')
c
c     ----- change orbsym to go from 0 to nsym-1 -----
c
      do 1 i=1,norbs
         orbsym(i)=orbsym(i)-1
    1 continue
c
c
      return
      end
