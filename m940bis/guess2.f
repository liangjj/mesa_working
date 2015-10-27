*deck @(#)guess2.f	1.1  11/30/90
      subroutine guess2(ptgues,nwks,ibuf,rbuf,lnbuf,
     $                  hguess,nwksg,diag,t1,eigvec,eigval,
     $                  ops,ntotal,repcor,prtflg,dagtyp,title)
c
c***begin prologue     guess2
c***date written       870813   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c
c***keywords           ci, guess
c***author             saxe, paul (lanl)
c***source             @(#)guess2.f	1.1   11/30/90
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
      character*(*) ops, title
      character*(*) prtflg
      character*80 htit
      integer ptgues(nwks)
      integer ibuf(2,lnbuf)
      integer nwks
      integer nwksg
      integer nnpg
      real*8 diag(nwks)
      real*8 rbuf(lnbuf)
      real*8 hguess(nwksg,nwksg)
      real*8 t1(nwks)
      real*8 eigvec(nwksg,nwksg)
      real*8 eigval(nwksg)
      real*8 eguess
      real*8 t
      real*8 repcor
c
      character*4 itoc
      integer intkey
      logical logkey, dagtyp
      dimension title(3)
c
      common /io/ inp,iout
c
c     ----- fill an array pointing to guess numbers, or zero if unused
c
      call iosys('read real '//title(3)//' from hamiltonian',
     1            nwks,t1,0,' ')
      call copy(t1,diag,nwks)
      call izero(ptgues,nwks)
c
      do 30 i=1,nwksg
         eguess=1.0d+30
         do 20 j=1,nwks
            if(t1(j).ge.eguess) go to 10
            refwlk=j
            eguess=t1(j)
 10         continue
 20      continue
         t1(refwlk)=1.0d+30
         ptgues(refwlk)=i
 30   continue
c
c     ----- read through the large hamiltonian and extract the small
c
      call rzero(hguess,nwksg*nwksg)
c
      npass=ntotal/lnbuf
      left = ntotal - npass*lnbuf
      call iosys('rewind all on hamiltonian',0,0,0,' ')
c
      do 50 pass=1,npass
         call iosys('read integer '//title(1)//' from hamiltonian'//
     $              ' without rewinding',2*lnbuf,ibuf,0,' ')
         call iosys('read integer '//title(1)//' from hamiltonian'//
     $              ' without rewinding',wptoin(lnbuf),rbuf,0,' ')
c
         do 40 n=1,lnbuf
            i=ptgues(ibuf(1,n))
            if (i.le.0) go to 40
            j=ptgues(ibuf(2,n))
            if (j.le.0) go to 40
            ii=max(i,j)
            jj=min(i,j)
            hguess(ii,jj)=hguess(ii,jj)+rbuf(n)
 40      continue
 50   continue
c
      if(left.ne.0) then
         call iosys('read integer '//title(1)//' from hamiltonian'//
     $              ' without rewinding',2*left,ibuf,0,' ')
         call iosys('read integer '//title(1)//' from hamiltonian'//
     $              ' without rewinding',wptoin(left),rbuf,0,' ')
         do 60 n=1,left
            i=ptgues(ibuf(1,n))
            if (i.le.0) go to 60
            j=ptgues(ibuf(2,n))
            if (j.le.0) go to 60
            ii=max(i,j)
            jj=min(i,j)
            hguess(ii,jj)=hguess(ii,jj)+rbuf(n)
 60      continue
      endif
c   
c     ----- add in the diagonals -----
c
      do 70 j=1,nwks
         i=ptgues(j)
         if (i.le.0) go to 70
         hguess(i,i)=diag(j)
 70   continue
c
      do 80 i=1,nwksg
         do 90 j=1,i
            hguess(j,i)=hguess(i,j)
 90      continue
 80   continue    
      if (logkey(ops,'print=ci=guess-matrix',.false.,' ')) then
        htit='guess hamiltonian'
        call prntrm(htit,hguess,nwksg,nwksg,nwksg,nwksg,iout)
      end if
c
c        ----- diagonalize -----
c
      call dsyev('v','l',nwksg,hguess,nwksg,eigval,t1,5*nwksg,info)
c
      if (info.ne.0) then
         call lnkerr('error in diagonalization')
      end if
c
      nroots=intkey(ops,'ci=nroots',1,' ')
      nroots=min(nwksg,nroots)
c
      if (prtflg.ne.'minimum') then
         write (iout,1) nwksg
      end if
c
      do 100 root=1,nwksg
c
c        ----- expand up the vector -----
c
         call rzero(t1,nwks)
c
         t=0.0d+00
         do 200 i=1,nwks
            if (ptgues(i).gt.0) then
               t1(i)=eigvec(ptgues(i),root)
               if (abs(t1(i)).gt.t) then
                  t=abs(t1(i))
                  ref=i
               end if
            end if
 200     continue
c
         if (prtflg.ne.'minimum') then
            write (iout,2) root,ref,diag(ref)+repcor,
     1                     eigval(root)+repcor, t1(ref)
         end if
c
         call iosys('write real "guess:'//itoc(root)//'" to guess',
     #        nwks,t1,0,' ')
 100  continue
      call iosys('rewind all on hamiltonian read-and-write',0,0,0,' ')
c
c
      return
 1    format (/,t15,'results from guess matrix of size',i5,/,
     1              ' guess   reference   diagonal energy  '
     2              'guess energy        c(0)')
 2    format (1x,i3,i10,4x,g20.9,g16.9,f8.4)
      end
