*deck @(#)guess2.f	5.1  11/6/94
      subroutine guess2(ptgues,nwks,ibuf,rbuf,lnbuf,hguess,nwksg,nnpg,
     $     d,t1,eigvec,eigval,ops,ntotal,repcor,prtflg)
c
c***begin prologue     guess2
c***date written       870813   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c
c***keywords           ci, guess
c***author             saxe, paul (lanl)
c***source             @(#)guess2.f	5.1   11/6/94
c
c***purpose            to form and diagonalize a portion of the h
c     matrix to obtain guesses at ci vectors.
c
c***description
c
c***references
c
c***routines called    (none)
c
c***end prologue       guess2
c
      implicit integer (a-z)
c
      character*(*) ops
      character*(*) prtflg
      integer ptgues(nwks)
      integer ibuf(2,lnbuf)
      integer nwks
      integer nwksg
      integer nnpg
      real*8 d(nwks)
      real*8 rbuf(lnbuf)
      real*8 hguess(nnpg)
      real*8 t1(nwks)
      real*8 eigvec(nwksg,nwksg)
      real*8 eigval(nwksg)
      real*8 eguess
      real*8 t
      real*8 repcor
c
      character*4 itoc
      integer intkey
      logical logkey
c
      common /io/ inp,iout
c
c     ----- fill an array pointing to guess numbers, or zero if unused
c
      call iosys('read real diagonals from hamiltonian',nwks,d,0,' ')
      call izero(ptgues,nwks)
c
      do 30 i=1,nwksg
         eguess=1.0d+30
         do 20 j=1,nwks
            if(d(j).ge.eguess) go to 10
            refwlk=j
            eguess=d(j)
 10         continue
 20      continue
         d(refwlk)=1.0d+30
         ptgues(refwlk)=i
 30   continue
c
c     ----- read through the large hamiltonian and extract the small
c
      call rzero(hguess,nnpg)
c
      npass=(ntotal-1)/lnbuf+1
      nleft=ntotal
      call iosys('rewind all on hamiltonian',0,0,0,' ')
c
      do 50 pass=1,npass
         call iosys('read integer buffers from hamiltonian'//
     $        ' without rewinding',2*lnbuf,ibuf,0,' ')
         call iosys('read integer buffers from hamiltonian'//
     $        ' without rewinding',wptoin(lnbuf),rbuf,0,' ')
c
         do 40 n=1,min(lnbuf,nleft)
            i=ptgues(ibuf(1,n))
            if (i.le.0) go to 40
            j=ptgues(ibuf(2,n))
            if (j.le.0) go to 40
            ii=max(i,j)
            jj=min(i,j)
            ij=ii*(ii-1)/2+jj
            hguess(ij)=hguess(ij)+rbuf(n)
 40      continue
         nleft=nleft-lnbuf
 50   continue
c
c     ----- add in the diagonals -----
c
      call iosys('read real diagonals from hamiltonian',nwks,d,0,' ')
      do 60 j=1,nwks
         i=ptgues(j)
         if (i.le.0) go to 60
         ii=i*(i+1)/2
         hguess(ii)=d(j)
 60   continue
c
      if (logkey(ops,'print=ci=guess-matrix',.false.,' ')) then
        call print(hguess,nnpg,nwksg,iout)
      end if
c
c        ----- diagonalize -----
c
      call rsp(nwksg,nwksg,nnpg,hguess,eigval,1,eigvec,d,t1,ierr)
c
      if (ierr.ne.0) then
         call lnkerr('error in rsp')
      end if
c
      nroots=intkey(ops,'ci=nroots',1,' ')
      nroots=min(nwksg,nroots)
      nguess=intkey(ops,'ci=nguess',nroots,' ')
c
      call iosys('read real diagonals from hamiltonian',nwks,d,0,' ')
      if (prtflg.ne.'minimum') then
         write (iout,70) nwksg
 70      format (/,t15,'results from guess matrix of size',i5,/,
     $        ' guess   reference   diagonal energy  guess energy',
     $        '   c(0)')
      end if
c
      do 90 root=1,nguess
c
c        ----- expand up the vector -----
c
         call rzero(t1,nwks)
c
         t=0.0d+00
         do 80 i=1,nwks
            if (ptgues(i).gt.0) then
               t1(i)=eigvec(ptgues(i),root)
               if (abs(t1(i)).gt.t) then
                  t=abs(t1(i))
                  ref=i
               end if
            end if
 80      continue
c
         if (prtflg.ne.'minimum') then
            write (iout,2) root,ref,d(ref)+repcor,eigval(root)+repcor,
     $           t1(ref)
 2          format (1x,i3,i10,4x,g20.9,g16.9,f8.4)
         end if
c
         call iosys('write real "guess:'//itoc(root)//'" to guess',
     #        nwks,t1,0,' ')
 90   continue
c
c
      return
      end
