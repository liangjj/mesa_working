*deck @(#)incore.f	1.3 7/30/91
      subroutine incore(lnbuf,ibuf,rbuf,nwks,nroots,eigvec,
     $                  eigval,t1,h,nnpwks,ops,rep,fzcore,
     $                  prtflg,calc,dagtyp,title)
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
      character*(*) ops, title
      character*(*) prtflg,calc
      integer ibuf(2,lnbuf)
      logical logkey, dagtyp
      character*80 hamtit
      real*8 rbuf(lnbuf)
      real*8 eigvec(nwks,nwks)
      real*8 eigval(nwks)
      real*8 t1(5*nwks)
      real*8 h(nwks,nwks)
      real*8 rep
      real*8 fzcore
      dimension title(3)   
c
      common /io/ inp,iout
      write(iout,1) nwks, nroots
c
c     ----- fill up the hamiltonian matrix -----
c
      call rzero(h,nwks*nwks)
      
c
      call iosys('read integer '//title(2)//' from hamiltonian',1,
     1           ntotal,0,' ')
c
      npass=ntotal/lnbuf
      left = ntotal - npass*lnbuf
      call iosys('rewind all on hamiltonian',0,0,0,' ')
c
      do 10 pass=1,npass
         call iosys('read integer '//title(1)//' from hamiltonian'//
     $              ' without rewinding',2*lnbuf,ibuf,0,' ')
         call iosys('read integer '//title(1)//' from hamiltonian'//
     $              ' without rewinding',wptoin(lnbuf),rbuf,0,' ')
c
         do 5 num=1,lnbuf
            i=ibuf(1,num)
            j=ibuf(2,num)
            h(i,j)=h(i,j)+rbuf(num)
 5       continue
 10   continue
c
      if(left.ne.0) then
         call iosys('read integer '//title(1)//' from hamiltonian'//
     $              ' without rewinding',2*left,ibuf,0,' ')
         call iosys('read integer '//title(1)//' from hamiltonian'//
     $              ' without rewinding',wptoin(left),rbuf,0,' ')
         do 20 num=1,left
            i=ibuf(1,num)
            j=ibuf(2,num)
            h(i,j)=h(i,j)+rbuf(num)
 20      continue
      endif
c
c     ----- add in the diagonals -----
c
      call iosys('read real '//title(3)//' from hamiltonian',
     1            nwks,t1,0,' ')
      do 30 i=1,nwks
         h(i,i) = h(i,i)+ t1(i)
 30   continue
      do 40 i=1,nwks
         do 50 j=1,i
            h(j,i)=h(i,j)
 50      continue
 40   continue   
c
      if (logkey(ops,'print=ci=hmatrix',.false.,' ')) then
         hamtit='hamiltonian matrix'
         call prntrm(hamtit,h,nwks,nwks,nwks,nwks,iout)
      end if
c
      if (logkey(ops,'print=ci=hmatrix=only',.false.,' ')) return
c
c        ----- diagonalize -----
c
      call dsyev('v','l',nwks,h,nwks,eigval,t1,5*nwks,info)
c
      if (info.ne.0) then
         call lnkerr('error in diagonalization')
      end if
c
c
      call iosys('read real '//title(3)//' from hamiltonian',
     1            nwks,t1,0,' ')
      call eig910(nwks,nroots,eigval,eigvec,t1,rep+fzcore,prtflg,
     #            calc)
c
c
      return
 1    format(/,15x,'direct diagonalization',/,/,5x,
     1             'size of matrix  = ',i3,/,5x,
     2             'number of roots = ',i3)
      end
