*deck besjyn
      subroutine besjyn(x,bj,dj,by,dy,n,npt,prn,nm)
c
c       ==========================================================
c       purpose: compute bessel functions jn(x) & yn(x) and
c                their derivatives
c       input :  x --- argument of jn(x) & yn(x)  ( x ò 0 )
c                n --- order of jn(x) & yn(x)
c       output:  bj(n) --- jn(x)
c                dj(n) --- jn'(x)
c                by(n) --- yn(x)
c                dy(n) --- yn'(x)
c                nm --- highest order computed
c       routines called:
c            (1) jy01b to calculate j0(x), j1(x), y0(x) & y1(x)
c            (2) msta1 and msta2 to calculate the starting 
c                point for backward recurrence
c       =========================================================
c
      implicit integer (a-z)
      real*8 x, bj, dj, by, dy
      real*8 bj0, dj0, bj1, dj1, by0, dy0, by1, dy1, bjk
      real*8 f, f0, f1, f2, cs, xinv 
      character*80 title
      logical prn
      dimension bj(npt,0:n), by(npt,0:n), dj(npt,0:n), dy(npt,0:n)
      dimension x(npt)
      common/io/ inp, iout
      do 100 pt=1,npt
         nm=n
         if (x(pt).lt.1.0d-100) then
             do 10 k=0,n
                bj(pt,k)=0.0d0
                dj(pt,k)=0.0d0
                by(pt,k)=-1.0d+300
                dy(pt,k)=1.0d+300
 10          continue   
             bj(pt,0)=1.0d0
             dj(pt,1)=0.5d0
         else
             call besjy01(x(pt),bj0,dj0,bj1,dj1,by0,dy0,by1,dy1)
             bj(pt,0)=bj0
             bj(pt,1)=bj1
             by(pt,0)=by0
             by(pt,1)=by1
             dj(pt,0)=dj0
             dj(pt,1)=dj1
             dy(pt,0)=dy0
             dy(pt,1)=dy1
             if (n.gt.1) then
                 xinv=1.d0/x(pt)
                 if (n.lt.int(0.9*x(pt))) then
                     twon=2
                     do 20 k=2,n
                         bjk=twon*xinv*bj1-bj0
                         bj(pt,k)=bjk
                         bj0=bj1
                         bj1=bjk
                         twon=twon+2
 20                  continue    
                 else
                     m=msta1(x(pt),200)
                     if (m.lt.n) then
                         nm=m
                     else
                         m=msta2(x(pt),n,15)
                     endif
                     f2=0.0d0
                     f1=1.0d-100
                     twon=m+m+2
                     do 30 k=m,0,-1
                        f=twon*xinv*f1-f2
                        if (k.le.nm) then
                            bj(pt,k)=f
                        endif
                        f2=f1
                        f1=f
                        twon=twon-2
 30                  continue   
                     if (abs(bj0).gt.abs(bj1)) then
                         cs=bj0/f
                     else
                         cs=bj1/f2
                     endif
                     do 40 k=0,nm
                         bj(pt,k)=cs*bj(pt,k)
 40                  continue    
                 endif
                 do 50 k=2,nm
                    dj(pt,k)=bj(pt,k-1)-k*xinv*bj(pt,k)
 50              continue   
                 f0=by(pt,0)
                 f1=by(pt,1)
                 twon=2
                 do 60 k=2,nm
                    f=twon*xinv*f1-f0
                    by(pt,k)=f
                    f0=f1
                    f1=f
                    twon=twon+2
 60              continue   
                 do 70 k=2,nm
                    dy(pt,k)=by(pt,k-1)-k*xinv*by(pt,k)
 70              continue   
             endif
         endif
 100  continue
      if(prn) then     
         title='regular bessel functions'
         call prntfm(title,bj,npt,nm+1,npt,nm+1,iout)
         title='first derivative of regular bessel functions'
         call prntfm(title,dj,npt,nm+1,npt,nm+1,iout)
         title='irregular bessel functions'
         call prntfm(title,by,npt,nm+1,npt,nm+1,iout)
         title='first derivative of irregular bessel functions'
         call prntfm(title,dy,npt,nm+1,npt,nm+1,iout)
      endif	 
      return
      end
