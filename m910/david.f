*deck @(#)david.f	1.3  7/30/91

      subroutine david(c,s,nwks,nguess,mxvc,root,nroots,dvdvec,dvdmat,
     $                 mxiter,nnpit,rep,fzcore,ops,ibuf,rbuf,lnbuf,
     $                 title,prtflg,calc)
c
c***begin prologue     david
c***date written       870801   (yymmdd)
c***revision date      910307   (yymmdd)
c
c    7 march    1991   rlm at lanl
c      adding argument to bliu to avoid mixed-mode problem in
c      the 'initialize' call.
c
c***keywords           ci, davidson's method
c***author             saxe, paul (lanl)
c***source             @(#)david.f	1.3   7/30/91
c
c***purpose            to find few lowest eigenvalues and eigenvectors of
c                      a real, symmetric matrix.
c***description
c
c***references
c
c***routines called    (none)
c
c***end prologue       david
c
      implicit integer (a-z)
c
      character*(*) ops, title
      character*(*) calc
      character*20 status
      character*4 itoc
      character*8 prtflg
c
      logical logkey, nodsk
c
      integer ibuf(2,lnbuf)
c
      real*8 rbuf(lnbuf)
      real*8 c(nwks,mxvc)
      real*8 s(nwks,mxvc)
      real*8 root(nroots)
      real*8 dvdvec(*)
      real*8 dvdmat(*)
      real*8 rep
      real*8 fzcore
      real*8 eguess
      real*8 eci
      real*8 cnverg
      real*8 thresh
      real*8 czero
      real*8 edav
      real*8 fpkey
      dimension title(3)
      data nodsk/.false./
c
      common /io/     inp,iout
c
c     ----- get the tolerances -----
c
      if(logkey(ops,'mcscf',.false.,' ')) then
       thresh=fpkey(ops,'mcscf=ci=tolerance',1.0d-07,' ')
       cnverg=fpkey(ops,'mcscf=ci=convergence',thresh,' ')
       thresh=fpkey(ops,'mcscf=ci=threshold',cnverg,' ')
       nattim=intkey(ops,'mcscf=ci=nroots-at-a-time',1,' ')
      else
       thresh=fpkey(ops,'ci=tolerance',1.0d-05,' ')
       cnverg=fpkey(ops,'ci=convergence',thresh,' ')
       thresh=fpkey(ops,'ci=threshold',cnverg,' ')
       nattim=intkey(ops,'ci=nroots-at-a-time',1,' ')
      endif
c
c
      if (prtflg.eq.'minimum') then
         status='noprint'
      else
         status='print'
      end if
c
      call iosys('read integer '//title(2)//' from hamiltonian',1,
     1            ntotal,0,' ')
      call iosys('read real '//title(3)//' from hamiltonian',
     $            nwks,c,0,' ')
      if(ntotal.le.lnbuf) then
         nodsk=.true.
      endif
c
      call bliu('initialize',status,c,thresh,nwks,mxiter,nroots,iout,
     #         nattim,rep+fzcore,cnverg,dvdvec)
c
      nvc=1
c      nguess=intkey(ops,'ci=nguess',nroots,' ')
c
      do 50 guess=1,nguess
         call iosys('read real "guess:'//itoc(guess)//'" from guess',
     $        nwks,c,0,' ')
c
c        ----- multiply the trial vector by h -----
c
         call hmult(ibuf,rbuf,lnbuf,title(1),ntotal,c,s,
     1              nwks,1,nodsk)
c
         call bliu('with vectors',0,c,s,nwks,mxiter,0,0,0,0,0,0)
 50   continue
c
c     ----- now iterate until we have extracted all the desired roots
c
 60   continue
         call bliu('solve',status,c,s,nwks,mxiter,0,0,0,dvdmat,root,
     #             dvdvec)
         if (status.eq.'converged') go to 90
         if(nguess.eq.nwks) go to 90
c
         nvc=0
 70      continue
            nvc=nvc+1
            call bliu('new trial',status,c(1,nvc),s(1,nvc),nwks,
     #                mxiter,0,0,0,0,root,dvdvec)
            if (status.ne.'done'.and.nvc.lt.mxvc) go to 70
            if (status.eq.'done'.and.nvc.eq.1) go to 60
            if (status.eq.'done') nvc=nvc-1
c
            call hmult(ibuf,rbuf,lnbuf,title(1),ntotal,c,s,
     1                 nwks,nvc,nodsk)
c
            do 80 ivc=1,nvc
               call bliu('with vectors',0,c(1,ivc),s(1,ivc),nwks,
     #                   mxiter,0,0,0,0,0,0)
 80         continue
            if (status.eq.'done') go to 60
            nvc=0
         go to 70
c
 90   continue
c
c     ----- recover each root's ci vector, find the most important
c           configuration, and write vector to rwf
c
      call iosys('read real diagonals from bliu',-1,s,0,' ')
c
      write(iout,*)' davidson diagonalization(david): prtflg ',
     # prtflg
      if (prtflg.ne.'minimum') then
         write (iout,100)
 100     format (///,' root   reference  guess energy    ci energy  ',
     #        '  davidson energy  c(0)')
      end if
c
      do 130 iroot=1,nroots
         call bliu('get vector',status,c,0,nwks,0,iroot,0,0,0,0,0)
         if (status.ne.'ok') go to 130
c..bhl 10/30/89
         if(calc.eq.'mcscf') then
         call iosys('write real "mc root '//itoc(iroot)//'" to rwf',
     #               nwks,c,0,' ')
         else
         call iosys('write real "ci root '//itoc(iroot)//'" to rwf',
     #               nwks,c,0,' ')
         end if
c..bhl 10/30/89
         czero=0.0d+00
         refwlk=0
         do 110 i=1,nwks
            if (abs(c(i,1)).gt.abs(czero)) then
               czero=c(i,1)
               refwlk=i
            end if
  110    continue
c
         eguess=s(refwlk,1)+rep+fzcore
         eci=root(iroot)+rep+fzcore
         edav=eci+(eci-eguess)*(1.0d+00-czero**2)
         if (iroot.eq.1) then
            call iosys('write real energy to rwf',1,
     $                  eci,0,' ')
         end if
c..bhl
      call iosys('write real "ci energy '
     $                        //itoc(iroot)//'" to rwf',1,eci,0,' ')
c..bhl
c
      if (prtflg.ne.'minimum') then
         write (iout,120) iroot,refwlk,eguess,eci,edav,czero
  120    format (1x,i3,i10,3g18.9,f8.4)
      end if
c
  130 continue
c
      write(iout,*)'  ci energies have been stored on rwf file'
      call bliu('finish',0,0,0,0,0, 0,0,0,0,0,0)
c
c
      return
      end
