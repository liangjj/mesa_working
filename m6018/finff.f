*deck @(#)finff.f	1.1 9/8/91
c***begin prologue     finff
c***date written       930502   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           finff, link 6018
c***author             schneider, barry (nsf)  
c***source             m6018
c***purpose            final free-free integrals
c***                   
c***references
c
c***routines called    iosys, util and mdutil
c***end prologue       finff
      subroutine finff (ovpp,ovpm,hpvhp,opp,opm,omm,hpp,nlm,nchan,
     1                  maxlm,dimc,words,prnt)
      implicit integer(a-z)
      real*8 omm
      complex*16 ovpp, ovpm, hpvhp, opp, opm, hpp
      character*3 itoc
      character*80 title
      logical prnt
      dimension ovpp(maxlm,nchan), ovpm(maxlm,nchan), hpvhp(maxlm,nchan)
      dimension opp(*), opm(*), omm(*), hpp(*)
      dimension nlm(dimc)
      common /io/ inp,iout
      locm=1
      do 10 ch=1,nchan
         call filff (ovpp(1,ch),ovpm(1,ch),hpvhp(1,ch),opp(locm),
     1               opm(locm),omm(locm),hpp(locm),maxlm,nlm(ch))
         if (prnt) then
             title='opp-'//itoc(ch)
             write (iout,1) (opp(ii), ii=locm,locm+nlm(ch)-1)
             title='opm-'//itoc(ch)
             write (iout,1) (opm(ii), ii=locm,locm+nlm(ch)-1)
             title='hpp-'//itoc(ch)
             write (iout,1) (hpp(ii), ii=locm,locm+nlm(ch)-1)
             title='omm-'//itoc(ch)
             write (iout,2) (omm(ii), ii=locm,locm+nlm(ch)-1)
         endif
         locm=locm+nlm(ch)
   10 continue
      words=locm-1                    
      return
    1 format(/,5x,'(',e15.8,1x,e15.8,')','(',e15.8,1x,e15.8,')')
    2 format(/,5x,1x,e15.8,1x,e15.8,1x,e15.8,1x,e15.8)
      end



