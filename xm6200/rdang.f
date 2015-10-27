*deck %W%  %G%
c***begin prologue     rdang
c***date written       930802   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           rdang, link m6200
c***author             schneider, barry (nsf)
c***source             m6200
c***purpose            read in angular information for separable
c***                   quadrature in theta and phi.
c
c***routines called
c***end prologue       rdang
      subroutine rdang(card,str,lmax,mmax,nthmax,nphmax,nang)
      implicit integer (a-z)
      character*(*) card, str
      character*30 thtyp, phtyp, chrkey
      real*8 pi, range 
      dimension range(2)
      data pi / 3.141592653589793d0  /
      common /io/ inp, iout      
      lmax=intkey(card,'maximum-l-value',2,' ')
      mmax=intkey(card,'maximum-m-value',2,' ')
      call iosys ('write integer "max l value '//str//'" to lamdat',
     1             1,lmax,0,' ')
      call iosys ('write integer "max m value '//str//'" to lamdat',
     1             1,mmax,0,' ')     
      thtyp=chrkey(card,'type-theta-quadrature','legendre',' ')
      phtyp=chrkey(card,'type-phi-quadrature','simpson',' ')
      call iosys('write character "theta quadrature type '//
     1            str//'" to lamdat',0,0,0,thtyp)
      call iosys('write character "phi quadrature type '//
     1            str//'" to lamdat',0,0,0,phtyp)
      ndef=(lmax+1)/2+1
      nthet=intkey(card,'theta-quadrature-order',ndef,' ')
      ndef=2*mmax+1
      nphi=intkey(card,'phi-quadrature-order',ndef,' ')
      nthmax=max(nthmax,nthet)
      nphmax=max(nphmax,nphi)
      call iosys('write integer "theta quadrature order '//
     1            str//'" to lamdat',1,nthet, 0,' ')  
      call iosys('write integer "phi quadrature order '//
     1            str//'" to lamdat',1,nphi,0,' ')
      nang=nthet*nphi
      write(iout,100) nthet, nphi 
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
  100 format(' theta quadrature order:',i3,2x,
     1       ' phi quadrature order:  ',i3)
  200 format(' theta range (',f15.8,',',f15.8,')')
  300 format(' phi range   (',f15.8,',',f15.8,')')      
      end

