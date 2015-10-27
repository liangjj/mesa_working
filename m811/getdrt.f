*deck @(#)getdrt.f	5.1  11/6/94
      subroutine getdrt(kadd,ladd,ijadd,ijgrp,iout,orbsym
     #,                 end,ngroup,nrefs,symorb,numij,nbf,norbs,
     #                  levfrm,nlevs,nrows,nsym)
c
c  read in drt information
c
      implicit integer (a-z)
c
      integer kadd(symorb),ladd(symorb),ijadd(numij),ijgrp(numij)
      integer iout(nbf),orbsym(norbs)
c
      call iosys('read integer kadd from rwf',-1,kadd,0,' ')
      call iosys('read integer ladd from rwf',-1,ladd,0,' ')
      call iosys('read integer ijadd from rwf',-1,ijadd,0,' ')
      call iosys('read integer ijgrp from rwf',-1,ijgrp,0,' ')
      call iosys('read integer iout from rwf',-1,iout,0,' ')
      call iosys('read integer orbsym from rwf',-1,orbsym,0,' ')
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
