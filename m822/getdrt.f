*deck @(#)getdrt.f	1.2  7/30/91
      subroutine getdrt(kadd,ladd,ijadd,ijgrp,iout,orbsym
     #,                 end,ngroup,nrefs,symorb,numij,nbf,norbs,
     #                  levfrm,nlevs,nrows,nsym,dsk)
c
c  read in drt information
c
      implicit integer (a-z)
      character*(*) dsk
c
      integer kadd(symorb),ladd(symorb),ijadd(numij),ijgrp(numij)
      integer iout(nbf),orbsym(norbs)
c
      call iosys('read integer kadd from '//dsk,-1,kadd,0,0)
      call iosys('read integer ladd from '//dsk,-1,ladd,0,0)
      call iosys('read integer ijadd from '//dsk,-1,ijadd,0,0)
      call iosys('read integer ijgrp from '//dsk,-1,ijgrp,0,0)
      call iosys('read integer iout from '//dsk,-1,iout,0,0)
      call iosys('read integer orbsym from '//dsk,-1,orbsym,0,0)
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
