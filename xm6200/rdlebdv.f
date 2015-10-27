*deck %W%  %G%
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
      subroutine rdlebdv(card,str,nleb,lmax,mmax,nthmax,nphmax,nang)
      implicit integer (a-z)
      character*(*) card, str
      real*8 pi, range 
      dimension range(2)
      data pi / 3.141592653589793d0  /
      common /io/ inp, iout      
c
      nleb=intkey(card,'lebedev-quadrature-order',17,' ')
      lmax=nleb/2
      mmax=nleb/2
      call iosys ('write integer "max l value '//str//'" to lamdat',
     1             1,lmax,0,' ')
      call iosys ('write integer "max m value '//str//'" to lamdat',
     1             1,mmax,0,' ')     
c
      call iosys ('write integer "lebedev quadrature order '//
     1             str//'" to lamdat',1,nleb,0,' ')
      if (nleb.eq.11) then
           nthet=50
           nphi=50
           nang=50
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
      range(1)=0.d0
      range(2)=180.d0
      call fparr(card,'theta-range',range,2,' ')
      range(1)=range(1)*pi/180.d0
      range(2)=range(2)*pi/180.d0 
      write(iout,200) range(1),range(2)
      call iosys('write real "theta range '//
     1            str//'" to lamdat',2,range,0,' ') 
      range(1)=0.d0
      range(2)=360.d0
      call fparr(card,'phi-range',range,2,' ')
      range(1)=range(1)*pi/180.d0
      range(2)=range(2)*pi/180.d0 
      write(iout,300) range(1),range(2)          
      call iosys('write real "phi range '//
     1            str//'" to lamdat',2,range,0,' ')                    
      return
  100 format(' generating lebedev quadrature of order = ',i3,/
     1       ' number of angular quadrature points = ',i4)    
  200 format(' theta range (',f15.8,',',f15.8,')')
  300 format(' phi range   (',f15.8,',',f15.8,')')      
      end

