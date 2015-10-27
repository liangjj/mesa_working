      subroutine ovtrans(hbbin,hbbout,nsym,nchan,nscat,nsch,
     1 tr,scr,nbfmax,nstate,iprint)
       implicit real*8 (a-h,o-z)
c
c  transforms bound-bound overlap  matrix to
c  orthogonal basis
c
      real*8 hbbin(nsym*(nsym+1)/2)
      real*8 hbbout(nbfmax,nbfmax,nstate)
      real*8 tr(nbfmax,nbfmax), scr(nbfmax,nbfmax,nchan)
      integer nscat(nchan),nsch(nbfmax,nchan)
c
      do 1 i=1,nbfmax
      do 1 j=1,nbfmax
      do 1 k=1,nstate
1     hbbout(i,j,k) = 0.0
c
c hbbin*tr
c
      do 100 ic=1,nchan
      nscic=nscat(ic)
      do 91 i=1,nsym
      do 91 j=1,nscic
      scr(i,j,ic) = 0.0
   91 continue
      do 80 k=1,nsym
      do 80 j=1,nscic
      do 80 i=1,nsym
      ii=max0(i,k)
      kk=i+k-ii
      ik=(ii-1)*ii/2+kk
80    scr(i,j,ic) = scr(i,j,ic) + hbbin(ik)*tr(k,nsch(j,ic))
100   continue
c
c tr*(hbbin*tr)
c
      do 200 ic=1,nchan
      nscic=nscat(ic)
      ist=ic*(ic+1)/2
      do 190 i=1,nscic
      do 190 j=1,nscic
      hbbout(i,j,ist)=ddot(nsym,tr(1,nsch(i,ic)),1,scr(1,j,ic),1)
190   continue
200   continue
      if(iprint.ne.0) then
      do 300 ic=1,nchan
      do 300 jc=1,ic
      write(6,107) ic,jc
107   format(//,' transformed bound-bound overlaps for channels:',2i4)
      ist=ic*(ic-1)/2+jc
      nsic=nscat(ic)
      nsjc=nscat(jc)
      do 300 isc=1,nsic
300   write(6,101) isc,(hbbout(isc,j,ist),j=1,nsjc)
101   format(1x,i3,3("(",f8.5,3x,f8.5,")",3x),/,
     &     (4x,3("(",f8.5,3x,f8.5,")",3x)))
      endif
      return
      end
