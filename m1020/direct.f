*deck @(#)direct.f	5.1  11/6/94
      subroutine direct(iwork,hessn,b0,work,nindep,nder,ops)
c
c***begin prologue     direct
c***date written       861213  (yymmdd)
c***revision date      930711  (yymmdd)
c   july 11, 1993      rlm at lanl
c      picking up hessian and right hand sides from tints as opposed to ints.
c***keywords           direct solution of cphf equations
c***author             saxe, paul (lanl)
c***source             @(#)direct.f	5.1   11/6/94
c***purpose            direct solution of cphf equations
c***description
c
c***references
c***routines called
c***end prologue       direct
c
      implicit integer (a-z)
c
      character*(*) ops
      real*8 hessn(nindep,nindep),b0(nindep,nder)
      real*8 work(nindep)
      integer iwork(nindep)
      logical logkey
c
      common /io/ inp,iout
c
c     ----- read in the hessian and right hand sides -----
c
      call iosys('read real "scf hessian" from tints',nindep**2,hessn,
     #            0,' ')
      call iosys('read real b0 from tints',nindep*nder,b0,0,' ')
c
c     ----- solve the first equation -----
c
      call sgefs(hessn,nindep,nindep,b0(1,1),1,error,work,iwork)
c
      write (iout,1) error
    1 format (' m1020:cphf solution by direct method',/,
     #        5x,'digits accuracy of result:',i3)
c
      if (error.lt.0) call lnkerr('problems directly solving cphf')
c
c     ----- and pass through the rest of the equations -----
c
      do 2 der=2,nder
         call sgefs(hessn,nindep,nindep,b0(1,der),2,error,work,iwork)
    2 continue
c
      call iosys('write real "independent cphf solutions" to rwf',
     $     nindep*nder,b0,0,' ')
c
c     ----- print the solutions if requested -----
c
      if (logkey(ops,'print=cphf=solutions',.false.,' ')) then
         write (iout,3)
    3    format (/,t6,'independent cphf solutions')
         call matout(b0,nindep,nder,nindep,nder,iout)
      end if
c
c
      return
      end
