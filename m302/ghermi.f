*deck ghermi.f	
      subroutine ghermi(h,wt,b,maxrts,mxpts,prnt)
c
c***purpose: to calculate gauss-hermite points and weights
c
c barry schneider              4 april 2001
c
      implicit integer (a-z)
      real*8 h, wt, b, endpts
      logical prnt
      character*80 title
      character*3 itoc
      dimension h(mxpts), wt(mxpts), b(maxrts), endpts(2)
      common /io/inp, iout
c
      count=0
      do 20 i=1,maxrts
         call gaussq('hermite',i,0.d0,0.d0,0,endpts,b,
     1                h(count+1),wt(count+1))
         count=count+i
 20   continue   
      if(prnt) then
         count=0
         do 30 i=1,maxrts
            title='points and weights for '//itoc(i)
     1            //' Gauss-Hermite quadrature'
            write(iout,1)
            do 40 j=1,i
               count=count+1
               write(iout,2) h(count), wt(count)
 40         continue
 30      continue
      endif   
      return
 1    format(/,6x,'     point     ','    weight     ')
 2    format(5x,2e15.8)
      end
