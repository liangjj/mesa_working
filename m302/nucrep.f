*deck @(#)nucrep.f	5.1  11/6/94
      subroutine nucrep(zan,c,nat,ncharge)
c***begin prologue     nucrep.f
c***date written       840820  
c***revision date      11/6/94      
c     25 january, 1994   rlm at lanl
c        modifying routine to recognize point charge background fields.
c        for solvent calculations. we want the interaction of the
c        point charges with the atomic nuclei, but not the
c        point charge-point charge interaction.
c
c***keywords           
c***author             saxe, paul(lanl)
c***source             @(#)nucrep.f	5.1   11/6/94
c***purpose            
c***description
c   to calculate the nuclear repulsion energy given the
c   charge of the nuclei, their locations and the number
c   of nuclei. 
c     
c***references
c
c***routines called
c
c***end prologue       nucrep.f
      implicit none
c     --- input variables -----
      integer nat,ncharge
c     --- input arrays (unmodified) ---
      real*8 zan(nat+ncharge),c(3,nat+ncharge)
c     --- input arrays (scratch) ---
c     --- output arrays ---
c     --- output variables ---
c     --- scratch arrays ---
c     --- local variables ---
      integer inp,iout
      integer i,j
      real*8 small,enuc,r
c
      common /io/     inp,iout
c
c     --- calculate the repulsion energy ---
      enuc=0.0d+00
      small=1.0d-6
      do 2 i=2,nat+ncharge
         do 1 j=1,i-1
            if(j.le.nat) then
               if(abs(zan(i)).gt.small.and.abs(zan(j)).gt.small) then
                  r=sqrt((c(1,i)-c(1,j))**2+(c(2,i)-c(2,j))**2+
     $                   (c(3,i)-c(3,j))**2)
                  if (r.lt.small) then
                     write (iout,3) i,j
    3                format (//,' $$$$$ nucrep: atoms',i3,' and',
     $                            i3,' are too close',//)
                     call lnkerr(' ')
                  end if
                  enuc=enuc+zan(i)*zan(j)/r
               endif
            endif
    1    continue
    2 continue
c
c     --- put the repulsion in constants ---
      call iosys('write real "nuclear repulsion energy" to rwf',
     $           1,enuc,0,' ')
c
c
      return
      end
