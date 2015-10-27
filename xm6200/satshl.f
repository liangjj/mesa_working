*deck %W%  %G%
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
      integer cskipb
      dimension nr(numshl), r(numshl+1)
      common /io/ inp, iout
      character*(*) type
      character*30 str
      character*800 card
      character*3 itoc, chrv
      logical logkey, nonsep,positn
      if (type.eq.'atom') then
          write(iout,100) icen      
          chrv=itoc(icen)
          len=cskipb(chrv,' ')
          str='atom-'//chrv(1:len)
          len=cskipb(str,' ')
      elseif (type.eq.'scattering') then
          write(iout,200)
          str='scattering center'
          len=cskipb(str,' ')
      else
          call lnkerr('error in center type')
      endif
      if (positn('$quadrature-'//str(1:len),card,inp)) then
         call cardin(card)
         nonsep=logkey(card,'angular-quadrature=lebedev',.false.,' ')
         if (.not.nonsep) then
             nleb=0
             write(iout,300)
             call iosys('write character "non separable quadrature" '//
     1                  'to lamdat',0,0,0,'no')
             call rdang(card,str,lmax,mmax,nthmax,nphmax,nang)
         else
            write(iout,400)
            call iosys('write character "non separable quadrature" '//
     1                 'to lamdat',0,0,0,'yes')
            call rdlebdv(card,str,nleb,lmax,mmax,nthmax,nphmax,nang)     
         endif
         ltop=max(ltop,lmax)
         mtop=max(mtop,mmax)
         call radquad(card,str,nr,r,numshl,igrid,iwt,nrmax,
     1                maxgrd,maxshl,ntrad,nang)
      else
         call lnkerr('cannot find the quadrature section')
      endif
c
c
  100 format(/,25x,' quadrature data for atom = ',i2)
  200 format(/,25x,' quadrature data for scattering center')
  300 format(' angular quadrature is separable in angles')
  400 format(' angular quadrature is non-separable (lebedev) in angles')
      return
      end







