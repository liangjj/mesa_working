*deck @(#)nucrep.f	5.1  11/6/94
      subroutine nucrep(zan,c,nat,grad)
c
c***purpose: to calculate the nuclear repulsion energy given the
c            charge of the nuclei, there locations and the number
c            of nuclei. the repulsion energy is stored on tape10.
c
c paul saxe               20 august 1984                          lanl
c
c 29 november 1985   pws @ lanl  modified for gradients
c
      implicit integer (a-z)
c
      real*8 zan(nat),c(3,nat),enuc,r,grad(3,nat)
c
      common /io/     inp,iout
c
c     ----- start timing -----
c
c
c     ----- calculate the repulsion energy -----
c
      enuc=0.0d+00
      call rzero(grad,3*nat)
      do 2 i=2,nat
         do 1 j=1,i-1
            r=sqrt((c(1,i)-c(1,j))**2+(c(2,i)-c(2,j))**2+
     #             (c(3,i)-c(3,j))**2)
            if (r.lt.1.0d-06) then
               write (iout,3) i,j
    3          format (//,' ##### nucrep: atoms',i3,' and',i3,' are',
     #                 ' too close',//)
               call lnkerr(' ')
            end if
            enuc=enuc+zan(i)*zan(j)/r
            do 9 coord=1,3
               grad(coord,i)=grad(coord,i)-(c(coord,i)-c(coord,j))*
     #                           zan(i)*zan(j)/r**3
               grad(coord,j)=grad(coord,j)+(c(coord,i)-c(coord,j))*
     #                           zan(i)*zan(j)/r**3
    9       continue
    1    continue
    2 continue
c
c     write (iout,4) enuc
c   4 format (/,' the nuclear repulsion energy is ',f15.9)
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
