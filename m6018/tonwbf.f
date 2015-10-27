*deck @(#)tonwbf.f	1.1 9/8/91
c***begin prologue     tonwbf
c***date written       930417   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           tonwbf, link 6018, kohn variational
c***author             schneider, barry (lanl)  
c***source             m6018
c***purpose            transform channel labelled bound-free integrals
c***                   to basis diagonalizing kinetic energy operator
c***
c***references
c***routines called    iosys, util and mdutil
c***end prologue       tonwbf
      subroutine tonwbf (ovbfi,hpvbi,ovbfo,hpvbo,keout,ovpb,hpb,
     1                   nbtot,orblst,nlm,nchan,maxlm,nmo,
     2                   dimmo,dimc)
      implicit integer(a-z)
      real*8 keout
      complex*16 ovbfi, hpvbi, ovbfo, hpvbo, ovpb, hpb
      dimension ovbfi(maxlm,nmo,nchan), hpvbi(maxlm,nmo,nchan)
      dimension ovbfo(*), hpvbo(*), ovpb(*), hpb(*)
      dimension keout(*), nbtot(dimc), orblst(dimmo,dimc), nlm(dimc)
      common /io/ inp,iout
      loch=1
      loco=1
      do 10 ch=1,nchan
c----------------------------------------------------------------------c
c            fill the temporary matrices ovpb and hpb in packed        c
c            form to prepare for transformation.                       c
c----------------------------------------------------------------------c
         call ovfil(ovbfi(1,1,ch),hpvbi(1,1,ch),ovpb,hpb,orblst,
     1              nlm(ch),nbtot(ch),nmo,maxlm,dimmo)
c----------------------------------------------------------------------c
c            transform and put back into ovbfo and hpvbo but in        c
c            packed form. the original dimensioning of these matrices  c
c            has enough space. remember only nbscat space has been     c
c            changed.                                                  c
c----------------------------------------------------------------------c
         call ecbc(ovbfo(loco),ovpb,keout(loch),nlm(ch),nbtot(ch),
     1             nbtot(ch))
         call ecbc(hpvbo(loco),hpb,keout(loch),nlm(ch),nbtot(ch),
     1             nbtot(ch))
         loch=loch+nbtot(ch)*nbtot(ch)
         loco=loco+nlm(ch)*nbtot(ch)
   10 continue
      return
      end



