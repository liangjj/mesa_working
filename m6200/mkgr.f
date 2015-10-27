*deck mkgr.f
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
      common /io/ inp, iout
      count=0
      if (nonsep) then
          nprd=nang
          do 10 i=1,nr
             do 20 j=1,nang
                z=rpt(i)*thpt(j)
                zz=rpt(i)*sthet(j)
                count=count+1
                grid(1,count)=zz*cphi(j)+cen(1)
                grid(2,count)=zz*sphi(j)+cen(2)
                grid(3,count)=z+cen(3)
   20        continue
   10     continue
      else
          nprd=nthet*nphi
          do 30 i=1,nr
             do 40 j=1,nthet
                z=rpt(i)*thpt(j)
                zz=rpt(i)*sthet(j)
                do 50 k=1,nphi
                   count=count+1
                   grid(1,count)=zz*cphi(k)+cen(1)
                   grid(2,count)=zz*sphi(k)+cen(2)
                   grid(3,count)=z+cen(3)
   50           continue
   40        continue
   30     continue
      endif
      do 60 icount=1,count
         rv=sqrt( grid(1,icount)*grid(1,icount) +
     1            grid(2,icount)*grid(2,icount) +
     2            grid(3,icount)*grid(3,icount) )
         costh=grid(3,icount)/rv
         if (abs(costh).ne.1.d0) then
             sinth=sqrt(1.d0-costh*costh)
             snphi=grid(2,icount)/(rv*sinth)
         else
             sinth=sqrt(1.d0-costh*costh)
             snphi=0.d0
         endif             
         spgrid(1,icount)=rv
         spgrid(2,icount)=costh
         spgrid(3,icount)=asin(snphi)
   60 continue      
      return
      end

