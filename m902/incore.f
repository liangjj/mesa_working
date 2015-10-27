*deck @(#)incore.f	5.1 11/6/94
      subroutine incore(lnbuf,ibuf,rbuf,nwks,nroots,eigvec,
     $                  eigval,t1,t2,h,nnpwks,ops,rep,fzcore,
     $                  prtflg,calc)
c
c***begin prologue     incore
c***date written       870730   (yymmdd)
c***revision date      910606   (yymmdd)
c
c   6 june     1991    rlm at lanl
c      passing nroots; useful for mcscf excited root optimizations.
c***keywords           ci diagonalization
c***author             saxe, paul (lanl)
c***source             @(#)incore.f	1.1   11/30/90
c
c***purpose            to brute-force diagonalize a ci hamiltonian
c
c***description
c
c***references
c
c***routines called    (none)
c
c***end prologue       incore
c
      implicit integer (a-z)
c
      character*(*) ops
      character*(*) prtflg,calc
      integer ibuf(2,lnbuf)
      logical logkey
      real*8 rbuf(lnbuf)
      real*8 eigvec(nwks,nwks)
      real*8 eigval(nwks)
      real*8 t1(nwks)
      real*8 t2(nwks)
      real*8 h(nnpwks)
      real*8 rep
      real*8 fzcore
c
      common /io/ inp,iout
c
c     ----- fill up the triangular hamiltonian matrix -----
c
      call rzero(h,nnpwks)
c
      call iosys('read integer "number of elements" from hamiltonian',
     $           1,ntotal,0,' ')
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
            ij=ibuf(1,i)*(ibuf(1,i)-1)/2+ibuf(2,i)
            h(ij)=h(ij)+rbuf(i)
 5       continue
         ntotal=ntotal-lnbuf
 10   continue
c
c     ----- add in the diagonals -----
c
      call iosys('read real diagonals from hamiltonian',nwks,t1,0,' ')
      do 20 i=1,nwks
         ii=i*(i+1)/2
         h(ii)=t1(i)
 20   continue
c
      if (logkey(ops,'print=ci=hmatrix',.false.,' ')) then
         call print(h,nnpwks,nwks,iout)
      end if
c
      if (logkey(ops,'print=ci=hmatrix=only',.false.,' ')) return
c
c        ----- diagonalize -----
c
      call rsp(nwks,nwks,nnpwks,h,eigval,1,eigvec,t1,t2,ierr)
c
      if (ierr.ne.0) then
         call lnkerr('error in rsp')
      end if
c
c
      call iosys('read real diagonals from hamiltonian',nwks,t1,0,' ')
      call eig902(nwks,nroots,eigval,eigvec,t1,rep+fzcore,prtflg,
     #            calc)
c
c
      return
      end
