*deck odeint.f 
c***begin prologue     odeint
c***date written       980924   (yymmdd)
c***revision date               (yymmdd)
c***keywords           odeint
c***                   
c***author             numerical recipe
c***source             
c***purpose            integrate a set of coupled differential equations
c***                   using a fifth order runge-kutta, variable stepsize
c***                   method.  
c
c***description        method is described in numerical recipes
c***references       
c
c***routines called   
c***end prologue       odeint
c
      subroutine odeint(ystart,x,yp,xp,x1,x2,dxsav,eps,
     1                  h1,hmin,t,n,nok,nbad,kmax,count,nstp)
c
c     t is a temporary variable allocated here which is used for
c     storage of the arrays y, dy, and yscal in the original code
c     as well as temporary storage for rkqs and rkck
c 
      implicit integer (a-z)
      real*8 ystart, x1, x2, yp, xp, t, dxsav, tiny
      real*8 eps, h1, hmin, h, hdid, hnext, x, xsav, xtest, dum
      dimension ystart(n), yp(n,*), xp(*), t(n,11)
      common/io/ inp, iout
      parameter ( tiny = 1.d-30 )    
      x=x1
      xtest=(x-x2)*(x2-x1)
      h=sign(h1,x2-x1)
      nok=0
      nbad=0
      count=0
      do 10 i=1,n
        t(i,1)=ystart(i)
 10   continue
      if (kmax.gt.0) then
          xsav=x-2.d0*dxsav
      endif
      do 20 stp=1,nstp
c      do while(xtest.ge.0.d0)
         call df(x,t(1,1),t(1,2),dum,dum,n)
         do 30 i=1,n
            t(i,3)=abs(t(i,1))+abs(h*t(i,2))+tiny
 30      continue
         if(kmax.gt.0)then
            if(abs(x-xsav).gt.abs(dxsav)) then
               if(count.lt.kmax-1)then
                  count=count+1
                  xp(count)=x
                  call copy(t,yp(1,count),n)
                  xsav=x
               endif
            endif
         endif
         if((x+h-x2)*(x+h-x1).gt.0.d0) then
             h=x2-x
         endif
         call rkqs(t(1,1),t(1,2),t(1,3),t(1,4),t(1,5),t(1,6),
     1             x,h,eps,hdid,hnext,n)
         if(hdid.eq.h)then
            nok=nok+1
         else
            nbad=nbad+1
         endif
         write(iout,*) x
         xtest=(x-x2)*(x2-x1)
         if(xtest.ge.0.d0) then
             call copy(t,ystart,n)
             if(kmax.ne.0)then
                count=count+1
                xp(count)=x
                call copy(t,yp(1,count),n)
             endif
             return
         endif
         if(abs(hnext).lt.hmin) then
c           write(iout,1) hnext
           call lnkerr('stepsize smaller than minimum in odeint')
            h=hnext
         endif
c      enddo
 20   continue
      call copy(t(1,1),ystart,n)
      write(iout,2) nstp, stp, count
      return
 1    format(/,1x,'next stepsize = ',e15.8)
 2    format(/,1x,'requested number of steps = ',i5,/,1x,
     1            'number of steps taken     = ',i5,/,1x,
     2            'number of steps stored    = ',i5 )
      end


