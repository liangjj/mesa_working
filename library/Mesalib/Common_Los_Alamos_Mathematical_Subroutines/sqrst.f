*deck sqrst
      subroutine sqrst (x, ldx, n, p, y, tol, b, rsd, k, 
     1                  jpvt, qraux, work)
      integer ldx, n, p, k, jpvt, info, j, kk
      real*8 x, y, tol, b, rsd, qraux, t
      dimension x(ldx,p), y(n), jpvt(p), qraux(p), work(p)
      dimension b(p), rsd(n)
      call izero(jpvt,p)
      call sqrdc(x,ldx,n,p,qraux,jpvt,work,1)
      k=0
      m=min(n,p)
      do 10 kk=1,m
         if(abs(x(kk,kk)).le.tol*abs(x(1,1))) go to 20
         k=kk
   10 continue
   20 continue   
      if(k.ne.0) then
         call sqrsl(x,ldx,n,k,qraux,y,rsd,rsd,b,rsd,rsd,110,info)
      endif
      do 30 j=1,p
         jpvt(j)=-jpvt(j)
         if(j.gt.k) then
            b(j)=0.d0
         endif
   30 continue
      do 40 j=1,p
         if(jpvt(j).gt.0) go to 40
         k=-jpvt(j)
   50    continue
         if(k.eq.j) go to 60
         t=b(j)
         b(j)=b(k)
         b(k)=t
         jpvt(k)=-jpvt(k)
         k=jpvt(k)
         go to 50
   60    continue
   40 continue
      return
      end
