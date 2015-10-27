*deck @(#)fmhpq.f	5.1  11/6/94
      subroutine fmhpq(lnbuf,ibuf,rbuf,hpq,nwksp,nwksq,ops)
c
c***begin prologue     fmhpq
c***date written       870916   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c
c***keywords           hpq
c***author             saxe, paul (lanl)
c***source             @(#)fmhpq.f	5.1   11/6/94
c
c***purpose            to form hpq from tape.
c
c***description
c
c***references
c
c***routines called    (none)
c
c***end prologue       fmhpq
c
      implicit integer (a-z)
c
      character*(*) ops
      integer ibuf(2,lnbuf)
      logical logkey
      real*8 rbuf(lnbuf)
      real*8 hpq(nwksp,nwksq)
c
      common /io/ inp,iout
c
c     ----- fill up the p-q block of the hamiltonian matrix -----
c
      call rzero(hpq,nwksp*nwksq)
c
      call iosys('read integer "number of elements" from hamiltonian',
     $     1,ntotal,0,' ')
c
      npass=(ntotal+lnbuf-1)/lnbuf
      call iosys('rewind all on hamiltonian',0,0,0,' ')
c
      do 10 pass=1,npass
         call iosys('read integer buffers from hamiltonian'//
     $        ' without rewinding',2*lnbuf,ibuf,0,' ')
         call iosys('read integer buffers from hamiltonian'//
     $        ' without rewinding',wptoin(lnbuf),rbuf,0,' ')
c
         do 5 i=1,min(lnbuf,ntotal)
            ii=ibuf(2,i)
            jj=ibuf(1,i)
            if (ii.gt.nwksp.or.jj.gt.nwksq) then
               call lnkerr('error in p,q section')
            end if
            hpq(ii,jj)=hpq(ii,jj)+rbuf(i)
 5       continue
         ntotal=ntotal-lnbuf
 10   continue
c
      if (logkey(ops,'print=ci=hpq',.false.,' ')) then
        write (iout,20)
 20     format(5x,'h(p,q)')
        call matout(hpq,nwksp,nwksq,nwksp,nwksq,iout)
      end if
c
c
      return
      end
