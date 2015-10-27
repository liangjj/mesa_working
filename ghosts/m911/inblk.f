*deck @(#)inblk.f	4.1  7/7/93
      subroutine inblk(aint,aone,nmax,norb2,itap49)
      character*6 caltyp,cityp,dertyp
      real*8 aint(nmax),aone(norb2),root(10)
      real*8 prtciv
      common /tapes/itap20,itape5,itape6,itape8,itap12,itap03,itap04
     *,             itape3,itap05,itap06
      common/dry2/prtciv,navail,nrused(2)
      data nwrite/0/
      save nwrite
      save nwrite
c
      do 20 i=1,norb2
         aone(i)=0.0d0
   20 continue
      return
c
c
      entry nxtblk
      do 10 i=1,nmax
         aint(i)=0.0d+00
   10 continue
      return
c
c
      entry putblk
      call iosys('write real "guga density matrix" to gden'//
     #' without rewinding',nmax,aint,0,' ')
c..bhl
c      call swrit(itap20,aint,intowp(nmax))
c..bhl
      return
c
c
      entry putone
c
      call iosys('write real "guga square 1pdm" to rwf without'//
     #' rewinding',norb2,aone,0,' ')
c..bhl
c      call srew(itap49)
c      call swrit(itap49,aone,intowp(norb2))
c..bhl
      return
      end
