*deck %W%  %G%
c***begin prologue     mkgr
c***date written       930802   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           mkgr, link m6200
c***author             schneider, barry (nsf)
c***source             m6200
c***purpose            calculate grid points in (x,y,z) and (r,cos,phi)
c***                   co-ordinates.
c***references         
c
c***routines called
c***end prologue       mkgr
      subroutine mkgr (grid,spgrid,rpt,thpt,sthet,sphi,cphi,cen,nr,
     1                   nthet,nphi,nang,nonsep)
      implicit integer (a-z)
      real*8 grid, spgrid, rpt, thpt, sthet, sphi, cphi, cen
      real*8 z, zz, rv, costh, sinth, snphi
      logical nonsep
      dimension grid(3,*), spgrid(3,*), rpt(nr), thpt(nthet)
      dimension sthet(nthet), sphi(nphi), cphi(nphi), cen(3)  
      real*8 zero,one,two,three,four,pi
      data zero/0.0d+00/
      data one/1.0d+00/,two/2.0d+00/,three/3.0d+00/,four/4.0d+00/
      save zero,one,two,three,four
c
      common /io/ inp, iout
      pi=four*atan(one)
      count=0
      if (nonsep) then
          do 10 i=1,nr
             do 20 j=1,nang
                z=rpt(i)*thpt(j)
                zz=rpt(i)*sthet(j)
                count=count+1
                grid(1,count)=zz*cphi(j)+cen(1)
                grid(2,count)=zz*sphi(j)+cen(2)
                grid(3,count)=z+cen(3)
                rv=sqrt( grid(1,count)*grid(1,count) +
     1                   grid(2,count)*grid(2,count) +
     2                   grid(3,count)*grid(3,count) )
                spgrid(1,count)=rv
                if(rv.eq.zero) then
                   spgrid(2,count)=zero
                   spgrid(3,count)=zero
                else
                   spgrid(2,count)=grid(3,count)/rv
                   if(grid(1,count).eq.zero) then
                      if(grid(2,count).ge.zero) then
                         spgrid(3,count)=pi/two
                      else
                         spgrid(3,count)=three*pi/two
                      endif
                   else
                      spgrid(3,count)=atan(grid(2,count)/grid(1,count))
                   endif
                endif
   20        continue
   10     continue
      else
          do 30 i=1,nr
             do 40 j=1,nthet
                z=rpt(i)*thpt(j)
                zz=rpt(i)*sthet(j)
                do 50 k=1,nphi
                   count=count+1
                   grid(1,count)=zz*cphi(k)+cen(1)
                   grid(2,count)=zz*sphi(k)+cen(2)
                   grid(3,count)=z+cen(3)
                   rv=sqrt( grid(1,count)*grid(1,count) +
     1                      grid(2,count)*grid(2,count) +
     2                      grid(3,count)*grid(3,count) )
                   spgrid(1,count)=rv
                   if(rv.eq.zero) then
                      spgrid(2,count)=zero
                      spgrid(3,count)=zero
                   else
                      spgrid(2,count)=grid(3,count)/rv
                      if(grid(1,count).eq.zero) then
                         if(grid(2,count).ge.zero) then
                            spgrid(3,count)=pi/two
                         else
                            spgrid(3,count)=three*pi/two
                         endif
                      else
                         spgrid(3,count)=
     $                      atan(grid(2,count)/grid(1,count))
                      endif
                   endif
   50           continue
   40        continue
   30     continue
      endif
      return
      end

