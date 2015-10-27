*deck @(#)modden.f	1.1 9/8/91
c***begin prologue     modden
c***date written       880423   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           modden, link 1106, kohn variational
c***author             schneider, barry (lanl)
c***source             m6008
c***purpose            modify optical potential denominator
c***              
c
c***references         
c
c***routines called    iosys, util and mdutil
c***end prologue       modden
      subroutine modden (hambb,energy,nbscat,nchan,nstri,mxb)
      implicit integer(a-z)
      real*8 energy, hambb
      dimension hambb(mxb,mxb,nstri), nbscat(nchan)
      do 10 i=1,nchan
         itri=i*(i+1)/2
         do 20 j=1,nbscat(i)
            hambb(j,j,itri)=hambb(j,j,itri)-energy
   20    continue
   10 continue
      return
      end
