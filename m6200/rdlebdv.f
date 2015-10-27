*deck rdlebdv.f
c***begin prologue     rdlebdv
c***date written       930802   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           rdlebdv, link m6200
c***author             schneider, barry (nsf)
c***source             m6200
c***purpose            read in lebedev quadrature information.
c
c***routines called
c***end prologue       rdlebdv
      subroutine rdlebdv(str,range,nleb,lmax,mmax,nthmax,nphmax,nang)
      implicit integer (a-z)
      character*30  str
      real*8 range 
      dimension range(4)
      common /io/ inp, iout      
      lmax=nleb/2
      mmax=nleb/2
      call iosys ('write integer "max l value '//str//'" to lamdat',
     1             1,lmax,0,' ')
      call iosys ('write integer "max m value '//str//'" to lamdat',
     1             1,mmax,0,' ')     
      call iosys ('write integer "lebedev quadrature order '//
     1             str//'" to lamdat',1,nleb,0,' ')
      if (nleb.eq.9) then
           nthet=38
           nphi=38
           nang=38
      elseif (nleb.eq.11) then
           nthet=50
           nphi=50
           nang=50
      elseif (nleb.eq.13) then
           nthet=74
           nphi=74
           nang=74
      elseif (nleb.eq.15) then
           nthet=86
           nphi=86
           nang=86
      elseif (nleb.eq.17) then
           nthet=110
           nphi=110
           nang=110
      elseif (nleb.eq.23) then
           nthet=194
           nphi=194
           nang=194
c      elseif (nleb.eq.29) then
c           nthet=302
c           nphi=302
c           nang=302
      else
           call lnkerr('lebedev quadrature not available for '//
     1                 'this order')
      endif
      write(iout,100) nleb, nang
      call iosys ('write integer "number lebedev points '//
     1             str//'" to lamdat',1,nang,0,' ')                         
      nthmax=max(nthmax,nthet)
      nphmax=max(nphmax,nphi)
      return
  100 format(' generating lebedev quadrature of order = ',i3,/
     1       ' number of angular quadrature points = ',i4)    
      end

