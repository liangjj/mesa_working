*deck @(#)mctdrt.f	5.1  11/6/94
      subroutine mctdrt(kadd,ladd,ijadd,ijgrp,bftorb,orbsym,
     $     end,ngroup,nrefs,symorb,numij,nbf,norbs,itape8,
     $     nijvir,ijww,klww,ijxx,klxx,levfrm,nlevs,nrows,
     $     nsym)
c
c***begin prologue
c***date written       871022   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c
c***keywords
c***author             lengsfield, byron (brl)
c***source             @(#)mctdrt.f	5.1   11/6/94
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
c  read in drt information
c
      implicit integer (a-z)
c
      integer kadd(symorb),ladd(symorb),ijadd(numij),ijgrp(numij)
      integer bftorb(nbf),orbsym(norbs)
      integer ijww(numij),klww(nijvir),ijxx(numij),klxx(nijvir)
c
c
      common /io/ inp,iout
c
c
      call iosys('read integer kadd from rwf',-1,kadd,0,' ')
      call iosys('read integer ladd from rwf',-1,ladd,0,' ')
      call iosys('read integer ijadd from rwf',-1,ijadd,0,' ')
      call iosys('read integer ijgrp from rwf',-1,ijgrp,0,' ')
cps      call iosys('read integer ningrp from rwf',-1,ningrp,0,' ')
      call iosys('read integer iout from rwf',-1,bftorb,0,' ')
      call iosys('read integer orbsym from rwf',-1,orbsym,0,' ')
      if (levfrm.gt.0) then
         call iosys('read integer ijxx from rwf',-1,ijxx,0,' ')
         call iosys('read integer klxx from rwf',-1,klxx,0,' ')
         call iosys('read integer ijww from rwf',-1,ijww,0,' ')
         call iosys('read integer klww from rwf',-1,klww,0,' ')
      end if
c
c
c     ----- change orbsym to go from 0 to nsym-1 -----
c
      do 1 i=1,norbs
         orbsym(i)=orbsym(i)-1
    1 continue
c
c      write (iout,2) ijxx
c    2 format (' ijxx',(t8,10i6))
c      write (iout,3) klxx
c    3 format (' klxx',(t8,10i6))
c      write (iout,4) ijww
c    4 format (' ijww',(t8,10i6))
c      write (iout,5) klww
c    5 format (' klww',(t8,10i6))
c
c
      return
      end
