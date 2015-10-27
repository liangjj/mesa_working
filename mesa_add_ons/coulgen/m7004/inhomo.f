*deck inhomo
      subroutine inhomo(rhs,psi0,psi1,sdrv,sumf,sumb,del,
     1                  fsumf,fsumb,ptmax,n)
c***begin prologue     inhomo
c***date written       920324   (yymmdd)
c***revision date               (yymmdd)
c***keywords
c***author             schneider, barry (nsf)
c***source             %W% %G% 
c***purpose            integrals from zero to x and from x to infinity
c***                                           i          i
c***                   are computed by a third order rule and used to
c***                   solve the integral equation for the driven
c***                   equation.
c***
c***description
c***            
c               
c               
c***references
c
c***routines called
c
c***end prologue       inhomo
c
      implicit integer (a-z)
      real*8 rhs, psi0, psi1, sdrv, del, wt, sumf, sumb, ff
      real*8 fsumf, fsumb
      dimension rhs(n), psi0(n), psi1(n), sdrv(n)
      dimension sumf(n), sumb(n)
      dimension wt(3,4)
      common /io/ inp, iout
      wt(1,1)=3.d0/8.d0
      wt(1,2)=19.d0/24.d0
      wt(1,3)=-5.d0/24.d0
      wt(1,4)=1.d0/24.d0
      wt(2,1)=1.d0/3.d0
      wt(2,2)=4.d0/3.d0
      wt(2,3)=wt(2,1)
      wt(2,4)=0.d0
      wt(3,1)=wt(1,1)
      wt(3,2)=3.d0*wt(3,1)
      wt(3,3)=wt(3,2)
      wt(3,4)=wt(3,1)  
      do 10 i=1,3
         do 20 j=1,4
            wt(i,j)=wt(i,j)*del
   20    continue
   10 continue
      sumf(1) = 0.d0
      do 30 i=1,n-3,3
         call vfill(sumf(i+1),sumf(i),3)
         cntj=0
         do 40 j=i,i+3
            ff=psi0(j)*rhs(j)
            cntj=cntj+1
            cntk=0
            do 50 k=i+1,i+3
               cntk=cntk+1
               sumf(k) = sumf(k) + wt(cntk,cntj)*ff
   50       continue
   40    continue
   30 continue
      sumb(n)=0.d0
      do 60 i=n,4,-3
         call vfill(sumb(i-3),sumb(i),3)
         cntj=0
         do 70 j=i,i-3,-1
            ff=psi1(j)*rhs(j)
            cntj=cntj+1
            cntk=0
            do 80 k=i-1,i-3,-1
               cntk=cntk+1
               sumb(k) = sumb(k) + wt(cntk,cntj)*ff
   80       continue 
   70    continue
   60 continue
      do 90 i=1,ptmax-1
         sumb(i)=sumb(i)-sumb(ptmax)
   90 continue
      sdrv(1)=psi0(1)*sumb(1)
      sdrv(ptmax)=psi1(ptmax)*sumf(ptmax)
      do 100 i=2,ptmax-1
         sdrv(i)=psi1(i)*sumf(i)+psi0(i)*sumb(i)
  100 continue
      fsumf=sumf(ptmax)
      fsumb=sumb(ptmax)
      ff=-1.d0/sumf(ptmax)
      do 200 i=1,ptmax
         sdrv(i)=sdrv(i)*ff
  200 continue
      return
      end















