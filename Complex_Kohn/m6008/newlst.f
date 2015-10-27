*deck @(#)newlst.f	1.1 9/8/91
c***begin prologue     newlst
c***date written       880423   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           newlst, link 1106, kohn variational
c***author             schneider, barry (lanl), rescigno, tom(llnl)  
c***source             m1106
c***purpose            index for bound matrix
c***description        forms index array for bound-bound and bound-free
c***                   matrix elements. variational orbitals first
c***                   for each channel. target orbitals at end by channel.
c***references         schneider and rescigno, physical review
c
c***routines called    iosys, util and mdutil
c***end prologue       newlst
      subroutine newlst (nbtot,nbscat,ntrgt,finlst,nchan,dimc,dimmo)
      implicit integer(a-z)
      dimension nbtot(dimc), finlst(dimmo,dimc), nbscat(dimc)
      dimension ntrgt(dimc)
      nscat=0
      do 100 i=1,nchan
  100 nscat=nscat+nbscat(i)
      orbcnt=0
      scatcn=nscat
      do 40 i=1,nchan
         do 30 j=1,nbscat(i)
            orbcnt=orbcnt+1
            finlst(j,i)=orbcnt
   30    continue
         if(ntrgt(i).ne.0) then     
            do 10 k=nbscat(i)+1,nbtot(i)
               scatcn=scatcn+1
               finlst(k,i)=scatcn
   10       continue
         endif
   40 continue
      return
      end
