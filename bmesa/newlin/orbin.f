      subroutine orbin (funct,pt,norb,lstorb,nsts,ldim,nptmx,n,ops)
c
      implicit integer(a-z)
      logical logkey
      common /io/ inp, iout
      real *8funct, pt
      character *3 itoc
      character *(*) ops
      dimension funct(nptmx,ldim,n), norb(nsts), lstorb(n,nsts)
      dimension pt(nptmx)
      nwtot=ldim*nptmx
*
*
c
c     ----- for each state read in the orbitals -----
c
      nmocnt=0
      do 30 is=1,nsts
      no=norb(is)
      if (no.eq.0) go to 30
      do 20 iorb=1,no
      offset=nwtot*(lstorb(iorb,is)-1)
      nmocnt=nmocnt+1
      call iosys ('read real "new s.c. mos" from mofile',nwtot,
     1            funct(1,1,nmocnt),offset,0)
   20 continue
   30 continue
      if (nmocnt.ne.n) then
      write (iout,70)
      stop 'orbs'
      endif
      if (logkey(ops,'print=lam=orbitals',.false.,' ')) then
      write (iout,50)
      do 40 i=1,nmocnt
      write (iout,60) i
      call matprt (funct(1,1,i),nptmx,ldim,nptmx,ldim,0,0,0,0,0,0,0)
   40 continue
      endif
      return
c
   50 format (/,5x,'orbitals')
   60 format (/,5x,'mo',1x,i4)
   70 format (//,5x,'inconsistency in no. of orbitals')
      end
