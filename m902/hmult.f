*deck @(#)hmult.f	5.1  11/6/94
         subroutine hmult(ibuf,rbuf,lnbuf,ntotal,c,s,nwks,nvecs)
c
c***begin prologue     hmult
c***date written       870801   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c
c***keywords           ci
c***author             saxe, paul (lanl)
c***source             @(#)hmult.f	5.1   11/6/94
c
c***purpose            multiply h by one or more c's using a
c                      diagonalization tape.
c
c***description
c
c***references
c
c***routines called    (none)
c
c***end prologue       hmult
c
      implicit integer (a-z)
c
      integer ibuf(2,lnbuf)
      real*8 rbuf(lnbuf)
      real*8 c(nwks,nvecs)
      real*8 s(nwks,nvecs)
      real*8 hij
c
c
      call rzero(s,nwks*nvecs)
c
      npass=(ntotal+lnbuf-1)/lnbuf
      nleft=ntotal
      call iosys('rewind all on hamiltonian',0,0,0,' ')
c
      do 30 pass=1,npass
         call iosys('read integer buffers from hamiltonian'//
     $        ' without rewinding',2*lnbuf,ibuf,0,' ')
         call iosys('read integer buffers from hamiltonian'//
     $        ' without rewinding',wptoin(lnbuf),rbuf,0,' ')
c
         do 20 n=1,min(lnbuf,nleft)
            i=ibuf(1,n)
            j=ibuf(2,n)
            hij=rbuf(n)
            do 10 vec=1,nvecs
               s(i,vec)=s(i,vec)+hij*c(j,vec)
               s(j,vec)=s(j,vec)+c(i,vec)*hij
 10         continue
 20      continue
         nleft=nleft-lnbuf
 30   continue
c
c
      return
      end
