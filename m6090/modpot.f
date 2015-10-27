*deck modpot
c***begin prologue     modpot
c***date written       920417   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           m6080, link 6080, model potential
c***author             schneider, barry (nsf)
c***source             m6080
c***purpose            to make mulichannel model potentials.
c***description        simple multichannel model potentials are set up
c***                   to test scattering calculation.
c***
c***
c***references         none
c
c***routines called
c***end prologue       modpot
      subroutine modpot(vmat,rgs,nptmx,charge,pot)
      implicit integer(a-z)
      real *8 vmat, rgs, vcpl, size, fpkey, charge, sfac
      character *800 card
      character *8 cpass
      character *16 pot, chrkey
      character *3 itoc
      dimension vmat(nptmx), rgs(1000)
      common /io/ inp, iout
      call posinp ('$pot',cpass)
      call cardin (card)
      pot=chrkey(card,'type','null',' ')
      size=fpkey(card,'cutoff',10.d0,' ')
      vcpl=fpkey(card,'barrier-size',1.d0,' ')
      sfac=fpkey(card,'shielding-factor',1.d0,' ')
      call rzero(vmat,nptmx)
      write(iout,100) pot
      if (pot.eq.'exponential') then
c     use exponential potential in all channels
      do 20 mp=1,nptmx
         vmat(mp)=-exp(-rgs(mp))
   20 continue
      elseif (pot.eq.'null') then
      do 60 mp=1,nptmx
         vmat(mp)=0.d+00
   60 continue
      elseif (pot.eq.'one') then
      do 70 mp=1,nptmx
         vmat(mp)=1.d+00
   70 continue
      elseif (pot.eq.'huck') then
      do 80 mp=1,nptmx
         vmat(mp)=vcpl
         if (rgs(mp).gt.size) vmat(mp)=0.d+00
   80 continue
      elseif (pot.eq.'shielded-coulomb') then
      do 200 mp=1,nptmx
         vmat(mp)=-exp(-sfac*rgs(mp))/rgs(mp)
  200 continue
      elseif(pot.eq.'coulomb') then
      do 90 mp=1,nptmx
         vmat(mp)=-charge/rgs(mp)
   90 continue
      else
      call lnkerr ('bad potential in modpot')
      endif
 100  format(//,10x,'potential type',a32)
      return
c
      end

