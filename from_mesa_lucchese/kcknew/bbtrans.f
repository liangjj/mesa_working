      subroutine bbtrans(hbbin,hbbout,nsym,nchan,nscat,nsch,
     1 tr,scr,nbfmax,nstate,iprint,hbbe,hbboutp,istat)
       implicit real*8 (a-h,o-z)
       character*8 istat,iex,ioft
c
c  reads direct hamiltonians from file made by mesa m950, then
c  transforms bound-bound hamiltonians  matrices to
c  orthogonal basis
c
c***note**** scr1 and hbbout are equivalenced in calling program
c
      real*8 hbbin(nsym**2),hbbe(*)
      real*8 hbbout(nbfmax,nbfmax,nstate),hbboutp(nbfmax,nbfmax,nstate)
      real*8 tr(nbfmax,nbfmax), scr(nbfmax,nbfmax,nstate)
      integer nscat(nchan),nsch(nbfmax,nchan)
      data iex/8hexchange/,ioft/8hoffshell/
c
c  read number of channels from mesa file bndmat as a check
c
      read(9) nroots
      if(nroots.ne.nchan) then
      write(6,10) nroots,nchan
  10  format(' trouble in kohnopt routine bbtrans',/,
     x '  number of channels from mesa .ne. number from ffbf etc.',
     x /,' values are: ',2i8)
      stop
      endif
c

      do 100 ic=1,nchan
      do 100 jc=1,ic
c
c read a hamiltonian for this channel pair
c
      read(9) mesaic,mesajc
      if(istat.eq.iex.or.istat.eq.ioft)then
      read(15) mesaic,mesajc
      call rdbinsqr(hbbe,nsym,nsym,15)
      endif
      call rdbinsqr(hbbin,nsym,nsym,9)
      if(iprint.ne.0) then
      write(6,*)' bound-bound direct hamiltonian'
      do 400 i=1,nsym
      k1=nsym*(i-1)+1
      k2=k1+nsym-1
 400     write(6,101)i,(hbbin(k),k=k1,k2)
      write(6,*)' bound-bound exchange hamiltonian'
      do 401 i=1,nsym
      k1=nsym*(i-1)+1
      k2=k1+nsym-1
 401     write(6,101)i,(hbbe(k),k=k1,k2)
      endif
      if(istat.eq.iex)then
      ij=0
      do 222 i=1,nsym
      do 222 j=1,nsym
      ij=ij+1
 222  hbbin(ij)=hbbin(ij)-.5*hbbe(ij)
      endif
c
c hbbin*tr
c
c
      nscjc=nscat(jc)
      ist=ic*(ic-1)/2+jc
      do 91 i=1,nsym
      do 91 j=1,nscjc
91    scr(i,j,ist)=0.
      do 80 k=1,nsym
      do 80 j=1,nscjc
      do 80 i=1,nsym
      ik=nsym*(k-1)+i
80    scr(i,j,ist) = scr(i,j,ist) + hbbin(ik)*tr(k,nsch(j,jc))
100   continue
c
c tr*(hbbin*tr)
c
      do 200 ic=1,nchan
      nscic=nscat(ic)
      do 200 jc=1,ic
      ist=ic*(ic-1)/2 + jc
      nscjc=nscat(jc)
      do 190 i=1,nscic
      do 190 j=1,nscjc
      hbbout(i,j,ist)=ddot(nsym,tr(1,nsch(i,ic)),1,scr(1,j,ist),1)
190   continue
200   continue
      if(istat.eq.ioft)then
      do 102 ic=1,nchan
      do 102 jc=1,ic
c
c hbbe*tr
c
c
      nscjc=nscat(jc)
      ist=ic*(ic-1)/2+jc
      do 92 i=1,nsym
      do 92 j=1,nscjc
 92      scr(i,j,ist)=0.
      do 81 k=1,nsym
      do 81 j=1,nscjc
      do 81 i=1,nsym
      ik=nsym*(k-1)+i
81    scr(i,j,ist) = scr(i,j,ist) + hbbe(ik)*tr(k,nsch(j,jc))
 102  continue
c
c tr*(hbbe*tr)
c
      do 201 ic=1,nchan
      nscic=nscat(ic)
      do 201 jc=1,ic
      ist=ic*(ic-1)/2 + jc
      nscjc=nscat(jc)
      do 191 i=1,nscic
      do 191 j=1,nscjc
      hbboutp(i,j,ist)=-.5*ddot(nsym,tr(1,nsch(i,ic)),1,scr(1,j,ist),1)
c      hbboutp(i,j,ist)=0.
 191  continue
201   continue
      endif
      if(iprint.ne.0) then
      do 300 ic=1,nchan
      do 300 jc=1,ic
      write(6,107) ic,jc
107   format(//,' transformed bound-bound ham. for channels:',2i4)
      ist=ic*(ic-1)/2+jc
      nsic=nscat(ic)
      nsjc=nscat(jc)
      do 300 isc=1,nsic
300   write(6,101) isc,(hbbout(isc,j,ist),j=1,nsjc)
101   format(1x,i3,6f12.5,/,(4x,6f12.5))
      if(istat.eq.ioft)then
      do 301 ic=1,nchan
      do 301 jc=1,ic
      write(6,108) ic,jc
 108  format(//,' transformed bound-bound exch ham. for channels:',2i4)
      ist=ic*(ic-1)/2+jc
      nsic=nscat(ic)
      nsjc=nscat(jc)
      do 301 isc=1,nsic
301   write(6,101) isc,(hbboutp(isc,j,ist),j=1,nsjc)
      endif
      endif
      return
      end
