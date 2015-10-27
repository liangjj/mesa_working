      subroutine modpot(varr,rgs,nptmx,charge,pot)
c***begin prologue     modpot
c***date written       861108   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           m1103, link 1103, model potential
c***author             schneider, barry (lanl)
c***source             m1103
c***purpose            to make mulichannel model potentials.
c***description        simple multichannel model potentials are set up
c***                   to test scattering calculation.
c***
c***
c***references         none
c
c***routines called
c***end prologue       modpot
      implicit integer(a-z)
      real *8 varr, rgs, vcpl, size, fpkey, charge, sfac
      character *800 card
      character *8 chrkey, pot, cpass
      character *3 itoc
      dimension varr(nptmx), rgs(1000)
      common /io/ inp, iout
      call posinp ('$inpot',cpass)
      call cardin (card)
      pot=chrkey(card,'type','exp',' ')
      size=fpkey(card,'cutoff',10.,' ')
      vcpl=fpkey(card,'barrier-size',1.,' ')
      sfac=fpkey(card,'shielding-factor',1.,' ')
      call rzero(varr,nptmx)
      if (pot.eq.'exp') then
c     use exponential potential in all channels
      do 20 mp=1,nptmx
         varr(mp)=-2.e+00*exp(-rgs(mp))
   20 continue
      elseif (pot.eq.'null') then
      do 60 mp=1,nptmx
         varr(mp)=0.e+00
   60 continue
      elseif (pot.eq.'one') then
      do 70 mp=1,nptmx
         varr(mp)=2.e+00
   70 continue
      elseif (pot.eq.'huck') then
      do 80 mp=1,nptmx
         varr(mp)=2.e+00*vcpl
         if (rgs(mp).gt.size) varr(mp)=0.d+00
   80 continue
      elseif (pot.eq.'shielded') then
      do 200 mp=1,nptmx
         varr(mp)=-2.*exp(-sfac*rgs(mp))/rgs(mp)
  200 continue
      elseif(pot.eq.'coulomb') then
      do 90 mp=1,nptmx
         varr(mp)=-2.e+00*charge/rgs(mp)
   90 continue
      else
      call lnkerr ('bad potential in modpot')
      endif
      return
c
      end
