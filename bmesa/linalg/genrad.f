c $Header$
*deck genrad.f
c***begin prologue     genrad
c***date written       921223   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           legndre, link 6203
c***author             schneider, barry (nsf)
c***source             m6203
c***purpose            generate radial solutions of free wave equation.
c***                   solve driven equation a la rescigno/orel for
c***                   smooth, cutoff irregular solution to driven wave
c***                   equation.
c***references         
c***routines called
c***end prologue       genrad
      subroutine genrad (greg,gireg,psln,gr1,gr2,energy,pt,wpt,rhsfn,
     1                   scr,wron,ns,nr,lval,nm,nl,nchan,lmax,
     2                   ldim,maxm,maxl,size)
      implicit integer(a-z)
      real*8 greg, gireg, gr1, gr2, energy, pt, wpt, rhsfn, scr, wron
      dimension greg(*), gireg(*), gr1(*), gr2(*), energy(nchan)
      dimension nr(ns), pt(*), wpt(*), nm(nchan), nl(maxm,nchan)
      dimension scr(*), wron(0:ldim), lval(maxl,maxm,nchan), rhsfn(*)
      character*2 itoc
      common /io/ inp, iout
c     read in the radial points and weights as one continuous vector
c                            then
c     check that size of array greg and gireg are big enough
      ptcnt=0
      asize=0
      do 10 is=1,ns
         pp=ptcnt+1
         call iosys ('read real "radial points shell-'//itoc(is)//
     1               '" from lamdat',nr(is),pt(pp),0,' ')
         call iosys ('read real "radial weights shell-'//itoc(is)//
     1               '" from lamdat',nr(is),wpt(pp),0,' ')
         do 20 chan=1,nchan
            do 30 mu=1,nm(chan)
               asize=asize+nl(mu,chan)*nr(is)
   30       continue
   20    continue
         ptcnt=ptcnt+nr(is)
   10 continue
      if (size.lt.asize) then
          call lnkerr('not enough space for greens functions')
      endif
      locrad=1
      do 40 chan=1,nchan
         do 50 mu=1,nm(chan)
c           calculate all l values from l = 0 to l = lmax using recursion
c           after these are computed pick the ones you want 
            call grnfn(gr1,gr2,energy(chan),pt,scr,ptcnt,lmax,
     1                 ldim,.false.)
            call gtherr(greg(locrad),gireg(locrad),gr1,gr2,
     1                  lval(1,mu,chan),ptcnt,nl(mu,chan),ldim)
c           solve driven free wave schroedinger equation
            call pertsl(psln(locrad),greg(locrad),gireg(locrad),pt,wpt,
     1                  rhsfn(locrad),ptcnt,nl(mu,chan),.false.)
            locrad=locrad+nl(mu,chan)*ptcnt
   50    continue
   40 continue
      return
      end
