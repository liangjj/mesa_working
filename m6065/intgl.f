*deck @(#)intgl.f
c***begin prologue     intgl
c***date written       xxxxxx   (yymmdd)
c***revision date      920409   (yymmdd)
c***keywords           kohn integrals
c***author             schneider, barry (nsf)
c***source             m6060
c***purpose            test integral
c***references
c
c***routines called    iosys, util and mdutil
c***end prologue       intgl
      subroutine intgl(sinnt,vlamda,grid,energy,npt,nolam)
      implicit integer (a-z)
      real *8 grid, energy, k, r, sn
      complex *16 sinnt, vlamda
      dimension vlamda(npt,nolam), sinnt(nolam), grid(4,npt)
c**********************************************************************c
c                the volume element is 4*pi*r*r                        c
c**********************************************************************c
      k=sqrt(energy)
      do 10 i=1,npt
         r=grid(1,i)*grid(1,i)+grid(2,i)*grid(2,i)+grid(3,i)*grid(3,i)
         r=sqrt(r)
         sn=sin(k*r)/r
         do 20 j=1,nolam
            sinnt(j)=sinnt(j)+sn*vlamda(i,j)
   20    continue
   10 continue     
      return
      end





