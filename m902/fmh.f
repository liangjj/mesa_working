*deck @(#)fmh.f	5.1  11/6/94
      subroutine fmh(lnbuf,ibuf,rbuf,h,nwks,t1,ops)
c
c***begin prologue     fmh
c***date written       870916   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c
c***keywords           h
c***author             saxe, paul (lanl)
c***source             @(#)fmh.f	5.1   11/6/94
c
c***purpose            to form h from tape.
c
c***description
c
c***references
c
c***routines called    (none)
c
c***end prologue       fmh
c
      implicit integer (a-z)
c
      character*(*) ops
      integer ibuf(2,lnbuf)
      logical logkey
      real*8 rbuf(lnbuf)
      real*8 h(nwks,nwks)
      real*8 t1(nwks)
c
      common /io/ inp,iout
c
c     ----- fill up the p-p block of the hamiltonian matrix -----
c
      call rzero(h,nwks**2)
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
            h(ibuf(1,i),ibuf(2,i))=h(ibuf(1,i),ibuf(2,i))+rbuf(i)
            h(ibuf(2,i),ibuf(1,i))=h(ibuf(2,i),ibuf(1,i))+rbuf(i)
 5       continue
         ntotal=ntotal-lnbuf
 10   continue
c
c     ----- add in the diagonals -----
c
      call iosys('read real diagonals from hamiltonian',nwks,t1,0,' ')
      do 15 i=1,nwks
         h(i,i)=t1(i)
 15   continue
c
      if (logkey(ops,'print=ci=hmatrix',.false.,' ')) then
        write (iout,20)
 20     format(5x,'h matrix')
        call matout(h,nwks,nwks,nwks,nwks,iout)
      end if
c
c
      return
      end
