c $Header$
*deck genleg.f
c***begin prologue     genleg
c***date written       921223   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           legndre, link 6203
c***author             schneider, barry (lanl)
c***source             m6203
c***purpose            generate real spherical harmonics from m = 0 to
c***                   mmax and l = m to lmax for each shell.
c***                   both m values and all l values are stored.      
c***references         
c***routines called
c***end prologue       genleg
      subroutine genleg (plm,phim,cphi,sphi,phi,wphi,theta,wthet,
     1                   dfct,ddfct,mpnt,lpnt,ns,nth,nph,lmax,
     2                   mmax,maxfac)
      implicit integer(a-z)
      dimension plm(*), phim(*), cphi(*), sphi(*), phi(*), theta(*)
      dimension wphi(*), wthet(*), mpnt(ns), lpnt(0:mmax,ns)
      dimension nth(ns), nph(ns)
      real*8 plm, phim, cphi, sphi, phi, theta, dfct, ddfct
      real*8 wphi, wthet
      character*2 itoc
      common /io/ inp, iout
c     mpnt is a pointer to the first word for m = 0 in each shell.
c     lpnt is a pointer to the first l value for a given m in each shell
      plmcnt=1
      phicnt=1
      phiwdc=1
      thewdc=1
      do 10 is=1,ns
         mpnt(is)=phicnt
         wds1=nph(is)+nph(is)
         call iosys ('read real "phi points shell-'//itoc(is)//
     1               '" from lamdat',nph(is),phi(phiwdc),0,' ')
         call iosys ('read real "phi weights shell-'//itoc(is)//
     1               '" from lamdat',nph(is),wphi(phiwdc),0,' ')
         call miscfn(phi(phiwdc),cphi,sphi,nph(is))
         call iosys ('read real "theta points shell-'//itoc(is)//
     1               '" from lamdat',nth(is),theta(thewdc),0,' ')
         call iosys ('read real "theta weights shell-'//itoc(is)//
     1               '" from lamdat',nth(is),wthet(thewdc),0,' ')
         do 20 m=0,mmax
            wds=(lmax+1-m)*nth(is)
            if (m.eq.0) then
                wds1=nph(is)
            endif
            lpnt(m,is)=plmcnt
            call scmphi(cphi,sphi,phim(phicnt),nph(is),m)
            call legend(plm(plmcnt),theta(thewdc),dfct,ddfct,nth(is),
     1                  lmax,m,maxfac)
            plmcnt=plmcnt+wds
            phicnt=phicnt+wds1
   20    continue
         phiwdc=phiwdc+nph(is) 
         thewdc=thewdc+nth(is)
   10 continue
      return
      end
