*deck satshl.f
c***begin prologue     satshl
c***date written       930802   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           satshl, link m6200
c***author             schneider, barry (nsf)
c***source             m6200
c***purpose            read in and write out atomic quadrature information
c***references         
c
c***routines called
c***end prologue       satshl
      subroutine satshl(type,icen,nr,r,numshl,igrid,iwt,nrmax,
     1                  nthmax,nphmax,maxgrd,maxshl,ntrad,ltop,
     2                  mtop,nleb,nang,nonsep)
      implicit real*8 (a-h,o-z)
      dimension nr(numshl), r(numshl+1), range(4)
      common /io/ inp, iout
      character*(*) type
      character*30 str, cpass
      character*800 card
      character*3 itoc, chrv
      logical logkey, nonsep
      data pi / 3.141592653589793d0  /
      if (type.eq.'atom') then
          write(iout,100) icen      
          chrv=itoc(icen)
          len=length(chrv)
          str='atom-'//chrv(1:len)
          len=length(str)
      elseif (type.eq.'scattering') then
          write(iout,200)
          str='scattering center'
          len=length(str)
      else
          call lnkerr('error in center type')
      endif
      call posinp('$quadrature-'//str(1:len),cpass)
      call cardin(card)
      nonsep=logkey(card,'angular-quadrature=lebedev',.false.,' ')
      if (.not.nonsep) then
          nleb=0
          write(iout,300)
          call iosys('write character "non separable quadrature" '//
     1               'to lamdat',0,0,0,'no')
          call rdang(card,str,lmax,mmax,nthmax,nphmax,nang)
      else
          write(iout,400)
          call iosys('write character "non separable quadrature" '//
     1               'to lamdat',0,0,0,'yes')
          nleb=intkey(card,'lebedev-quadrature-order',17,' ')
          range(1)=0.d0
          range(2)=180.d0
          call fparr(card,'theta-range',range,2,' ')
          range(1)=range(1)*pi/180.d0
          range(2)=range(2)*pi/180.d0 
          range(3)=0.d0
          range(4)=360.d0
          call fparr(card,'phi-range',range(3),2,' ')
          range(3)=range(3)*pi/180.d0
          range(4)=range(4)*pi/180.d0 
          write(iout,500) range(1),range(2)
          call iosys('write real "theta range '//
     1                str//'" to lamdat',2,range,0,' ') 
          write(iout,600) range(3),range(4)          
          call iosys('write real "phi range '//
     1                str//'" to lamdat',2,range(3),0,' ')  
          call rdlebdv(str,range,nleb,lmax,mmax,nthmax,nphmax,nang)     
      endif
      ltop=max(ltop,lmax)
      mtop=max(mtop,mmax)
      call radquad(card,str,nr,r,numshl,igrid,iwt,nrmax,
     1             maxgrd,maxshl,ntrad,nang)
  100 format(/,25x,' quadrature data for atom = ',i2)
  200 format(/,25x,' quadrature data for scattering center')
  300 format(' angular quadrature is separable in angles')
  400 format(' angular quadrature is non-separable (lebedev) in angles')
  500 format(' theta range (',f15.8,',',f15.8,')')
  600 format(' phi range   (',f15.8,',',f15.8,')')
      return
      end







