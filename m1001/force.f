*deck @(#)force.f	1.1  11/30/90
      subroutine force(ua,gama,bva,alagt,alag,ta,
     1                 ub,gamb,bvb,blagt,blag,tb,
     2  xlag, noc,nob, ncsf, xforc, prtall)
      implicit real*8 (a-h,o-z)
c
      common /io/ inp,iout
c
      dimension ua(nob,2),gama(2),bva(2),alagt(nob,2),
     1          alag(nob,2),ta(nob,2)
      dimension ub(nob,2),gamb(2),bvb(2),blagt(nob,2),
     1          blag(nob,2),tb(nob,2)
      dimension xlag(nob,2)
c
      logical prtall
c
c    alagt(j,i)*ub(j,i)
c    blag(j,i)*ta(j,i)
c
      xx=0.d0
      yy=0.d0
c
c
      do 10 i=1,noc
         do 20 j=1,nob
            xx=xx+alagt(j,i)*(ub(j,i)+tb(j,i))
            yy=yy+blag(j,i)*ta(j,i)
   20    continue
   10 continue
c
c    2*<c(i)|b(j)>
c
      zz=0.d0
      do 30 k=1,ncsf
         zz=zz+bva(k)*gamb(k)
   30 continue
c
c     scale and adjust phase
c
      zz=-2.d0*zz
c
c     lag(i,j)*(ti(j,k)*tj(i,k)+ti(k,j)*tj(k,j))
c
c     note lagrangian from mcscf needs to be doubled
c
      tt=0.d0
      do 40 i=1,noc
         do 50 j=1,noc
            x1=0.d0
            do 60 k=1,nob
               x1=x1+ta(k,i)*(tb(k,j)+tb(j,k))
   60       continue
            tt=tt+xlag(i,j)*x1
   50    continue
   40 continue
c
      tt=2.d0*tt
c
      if(prtall) then
         write(iout,70) xx,yy,zz,tt
  70     format(/,'  hessian contributions ',/,
     1            '   alagt*ub       ',f12.8,/,
     2            '   blag *ta       ',f12.8,/,
     3            '   2*<cb|ba>      ',f12.8,/,
     4            '   lag*(tb+tbt)ta ',f12.8)
      end if
c
      xforc=xx+yy+zz+tt
c
      return
      end
