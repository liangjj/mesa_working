*deck @(#)trnsf.f	1.1 9/9/91
      subroutine trnsf(grid,save,sum,suma,sumb,sumc,igrid,isave,
     1                 itran)
      parameter (nocen=10)
      implicit real *8 (a-h,o-z)
      dimension grid(4,igrid), save(isave,4)
      dimension ajacob(3,3),rq(3)
      logical itran, iplot
      character *16 ftitl
      common /centsi/ncent,iplot
      common /centsr/cc(3),a(3,nocen),eta(nocen),
     1              znuc(nocen),amu(nocen), p(nocen)
      common/io/ inp,iout
      data pi/3.14159265358d+00/
      iwrit=4*igrid
      ftitl='"untrns grid"'
      call iosys ('write integer "no. grid pts" to grid',1,igrid,0,' ')
      call iosys ('write real '//ftitl//' to grid',iwrit,grid,0,' ')
      sum = 0.0d+00
      suma = 0.0d+00
      sumb = 0.0d+00
      sumc = 0.0d+00
      do 40 ii=1,igrid
      xq=grid(1,ii)
      yq=grid(2,ii)
      zq=grid(3,ii)
      wnode=grid(4,ii)
      sum = sum+wnode
c
c
c  example is a yukawa potential at location a(1 to 3,i))
c
c
c
c  evaluate integrand in untransformed coordinates
c
      func = 0.d+00
      do 401 i=1,ncent
      dist =sqrt((xq-a(1,i))**2 + (yq-a(2,i))**2 + (zq-a(3,i))**2)
  401 func = func + exp(-eta(i)*dist)/dist
      suma = suma + func*wnode
      if (iplot) then
          if(ii.le.isave) then
             save(ii,2) = xq
             save(ii,4) = zq
          endif
      endif
c
c  evaluate the transformed coordinates themselves
c
          if (.not.itran) then
      xx = 0.d+00
      yy = 0.d+00
      zz = 0.d+00
      sump = 0.0d+00
      sump2=0.d+00
      do 402 i=1,ncent
      dist =sqrt((xq-a(1,i))**2 + (yq-a(2,i))**2 + (zq-a(3,i))**2)
      p(i) = exp(-amu(i)*dist**2)
      xx = xx+  (a(1,i)-xq)*p(i)**2
      yy = yy +  (a(2,i)-yq)*p(i)**2
      zz = zz +  (a(3,i)-zq)*p(i)**2
      sump = sump + p(i)
      sump2 = sump2 + p(i)**2
  402 continue
      if(sump.le.1.d-30) then
      xx = xq
      yy = yq
      zz = zq
      det = 1.d+00
      else
      xx = xq + xx/sump
      yy = yq + yy/sump
      zz = zq + zz/sump
c
c  evaluate jacobian of transformation
c
      rq(1) = xq
      rq(2) = yq
      rq(3) = zq
c
      do 420 n=1,3
      do 419 m=1,n
c
      sum1=0.d+00
      sum2=0.d+00
      sum3=0.d+00
      do 410 i = 1,ncent
      tmu = 2.d+00*amu(i)
      sum1=sum1+tmu*(rq(n)-a(n,i))*(rq(m)-a(m,i))*2.d+00*p(i)**2
      sum2=sum2+(rq(n)-a(n,i))*p(i)**2
      sum3=sum3+tmu*(rq(m)-a(m,i))*p(i)
 410   continue
c

      ajacob(n,m)=sum1/sump - (sum2/sump)*sum3/sump
      ajacob(m,n) = ajacob(n,m)
 419   continue
      ajacob(n,n)=ajacob(n,n)+1.d+00-sump2/sump
 420   continue
c
c evaluate determinant of jacobian
c
      det=    ajacob(1,1)*ajacob(2,2)*ajacob(3,3)
      det=det+ajacob(1,2)*ajacob(2,3)*ajacob(3,1)
      det=det+ajacob(2,1)*ajacob(3,2)*ajacob(1,3)
      det=det-ajacob(1,3)*ajacob(2,2)*ajacob(3,1)
      det=det-ajacob(1,2)*ajacob(2,1)*ajacob(3,3)
      det=det-ajacob(2,3)*ajacob(3,2)*ajacob(1,1)
c
c
c calculate integrand in transformed coordinates
c
      endif
      if (iplot) then
          if(ii.le.isave) then
              save(ii,1) = xx
              save(ii,3) = zz
          endif
      endif
      funca = 0.d0
      do 404 i=1,ncent
      dista=sqrt((xx-a(1,i))**2 + (yy-a(2,i))**2 + (zz-a(3,i))**2)
      funca = funca +  exp(-eta(i)*dista)/dista
404   continue
c
      wnode=abs(det)*wnode
      sumb=sumb+funca*wnode
      sumc = sumc + wnode
      endif
c
c load transformed points back into grid
c
      grid(1,ii)=xx
      grid(2,ii)=yy
      grid(3,ii)=zz
      grid(4,ii)=wnode
40    continue
c
c write transformed grid to disk
c
      ftitl='"trns grid"'
      call iosys ('write real '//ftitl//' to grid',iwrit,grid,0,' ')
      return
      end
