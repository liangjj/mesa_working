*deck  @(#)nucrep.f	5.1 11/6/94
      subroutine nucrep(zan,c,nat)
c
c***purpose: to calculate the nuclear repulsion energy given the
c            charge of the nuclei, there locations and the number
c            of nuclei. the repulsion energy is stored on tape10.
c
c paul saxe               20 august 1984                          lanl
c
      implicit integer (a-z)
c
      real*8 zan(nat),c(3,nat),small,enuc,r
c
      common /io/     inp,iout
c
c     ----- start timing -----
c
c
c     ----- calculate the repulsion energy -----
c
      enuc=0.0d+00
      small=1.d-6
      do 2 i=2,nat
         do 1 j=1,i-1
c..bhl 7/21/89 llnl
         if(abs(zan(i)).gt.small.and.abs(zan(j)).gt.small) then
            r=sqrt((c(1,i)-c(1,j))**2+(c(2,i)-c(2,j))**2+
     #             (c(3,i)-c(3,j))**2)
            if (r.lt.1.0d-06) then
               write (iout,3) i,j
    3          format (//,' ##### nucrep: atoms',i3,' and',i3,' are',
     #                 ' too close',//)
               call lnkerr(' ')
            end if
            enuc=enuc+zan(i)*zan(j)/r
          end if
c..bhl 7/21/89
    1    continue
    2 continue
c
c
c     ----- put the repulsion in constants -----
c
      call iosys('write real "nuclear repulsion energy" to rwf',
     $     1,enuc,0,' ')
c
c     ----- stop timing -----
c
c
c
      return
      end
