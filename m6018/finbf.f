*deck @(#)finbf.f	1.1 9/8/91
c***begin prologue     finbf
c***date written       930502   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           finbf, link 6018
c***author             schneider, barry (nsf)  
c***source             m6018
c***purpose            final bound-free integrals
c***                   
c***references
c
c***routines called    iosys, util and mdutil
c***end prologue       finbf
      subroutine finbf (ovbf,hpvb,opb,omb,hpb,hmb,nlm,nbtot,nchan,
     1                  dimc,words,prnt)
      implicit integer(a-z)
      real*8 omb, hmb
      complex*16 ovbf, hpvb, opb, hpb
      character*3 itoc
      character*80 title
      logical prnt
      dimension ovbf(*), hpvb(*), opb(*), omb(*), hpb(*), hmb(*)
      dimension nlm(dimc), nbtot(dimc)
      common /io/ inp,iout
      locm=1
      do 10 ch=1,nchan
         call filbf (ovbf(locm),hpvb(locm),opb(locm),omb(locm),
     1               hpb(locm),hmb(locm),nlm(ch),nbtot(ch))
         locm=locm+nlm(ch)*nbtot(ch)
         if (prnt) then
             title='opb-'//itoc(ch)
             call prntcm(title,opb(locm),nlm(ch),nbtot(ch),nlm(ch),
     1                   nbtot(ch),iout)
             title='hpb-'//itoc(ch)
             call prntcm(title,hpb(locm),nlm(ch),nbtot(ch),nlm(ch),
     1                   nbtot(ch),iout)
             title='omb-'//itoc(ch)
             call prntrm(title,omb(locm),nlm(ch),nbtot(ch),nlm(ch),
     1                   nbtot(ch),iout)
             title='hmb-'//itoc(ch)
             call prntrm(title,hmb(locm),nlm(ch),nbtot(ch),nlm(ch),
     1                   nbtot(ch),iout)
         endif
   10 continue
      words=locm-1                    
      return
      end



