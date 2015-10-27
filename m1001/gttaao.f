*deck @(#)gttaao.f	1.2  7/30/91
      subroutine gttaao(ta,scr,cm,nob,ldf)
c
      implicit integer(a-z)
      real*8 ta(*),scr(*),cm(*)
      common /io/ inp,iout
c
      nobnob=nob*nob
c
      call iosys('read real mo_der_overlap from dints'
     1           //' without rewinding',ldf*nobnob,ta,0,' ')
c..bhl
c        write(iout,*)' gttaao: mo orbitals '
c        call matout(cm,nob,nob,nob,nob,iout)
c..bhl
      ix=1
      do 10 i=1,ldf
c..bhl
c       if(i.eq.1) then
c        write(iout,*)' mo overlap ints ndf=1 '
c        call matout(ta,nob,nob,nob,nob,iout)
c       endif
c..bhl
         call fixta(ta(ix),nob)
c..bhl
c       if(i.eq.1) then
c          write(iout,*)' ta-mo  ndf=1 '
c          call matout(ta,nob,nob,nob,nob,iout)
c       endif
c..bhl
         call ebc(scr,cm,ta(ix),nob,nob,nob)
         call scopy(nobnob,scr,1,ta(ix),1)
c..bhl
c       if(i.eq.1) then
c          write(iout,*)' ta-ao ndf=1 '
c          call matout(ta,nob,nob,nob,nob,iout)
c       endif
c..bhl
        ix=ix+nobnob
  10  continue
c
      return
      end
