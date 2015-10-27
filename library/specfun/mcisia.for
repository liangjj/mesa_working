        program mcisia
c
c       ========================================================
c       purpose: this program computes the cosine and sine 
c                integrals using subroutine cisia
c       input :  x  --- argument of ci(x) and si(x)
c       output:  ci --- ci(x)
c                si --- si(x)
c       example:
c                    x         ci(x)          si(x)
c                 ------------------------------------
c                   0.0     - ì             .00000000
c                   5.0     -.19002975     1.54993124
c                  10.0     -.04545643     1.65834759
c                  20.0      .04441982     1.54824170
c                  30.0     -.03303242     1.56675654
c                  40.0      .01902001     1.58698512
c       ========================================================
c
        double precision ci,si,x
        write(*,*)'please enter x '
        read(*,*) x
        write(*,*)'   x         ci(x)          si(x)'
        write(*,*)'------------------------------------'
        call cisia(x,ci,si)
        if (x.ne.0.0d0) write(*,10)x,ci,si
        if (x.eq.0.0d0) write(*,20)
10      format(1x,f5.1,2f15.8)
20      format(3x,' .0',4x,' - ì',13x,'.00000000')
        end



        subroutine cisia(x,ci,si)
c
c       =============================================
c       purpose: compute cosine and sine integrals
c                si(x) and ci(x)  ( x ò 0 )
c       input :  x  --- argument of ci(x) and si(x)
c       output:  ci --- ci(x)
c                si --- si(x)
c       =============================================
c
        implicit double precision (a-h,o-z)
        dimension bj(101)
        p2=1.570796326794897d0
        el=.5772156649015329d0
        eps=1.0d-15
        x2=x*x
        if (x.eq.0.0d0) then
           ci=-1.0d+300
           si=0.0d0
        else if (x.le.16.0d0) then
           xr=-.25d0*x2
           ci=el+dlog(x)+xr
           do 10 k=2,40
              xr=-.5d0*xr*(k-1)/(k*k*(2*k-1))*x2
              ci=ci+xr
              if (dabs(xr).lt.dabs(ci)*eps) go to 15
10         continue
15         xr=x
           si=x
           do 20 k=1,40
              xr=-.5d0*xr*(2*k-1)/k/(4*k*k+4*k+1)*x2
              si=si+xr
              if (dabs(xr).lt.dabs(si)*eps) return
20         continue
        else if (x.le.32.0d0) then
           m=int(47.2+.82*x)
           xa1=0.0d0
           xa0=1.0d-100
           do 25 k=m,1,-1
              xa=4.0d0*k*xa0/x-xa1
              bj(k)=xa
              xa1=xa0
25            xa0=xa
           xs=bj(1)
           do 30 k=3,m,2
30            xs=xs+2.0d0*bj(k)
           bj(1)=bj(1)/xs
           do 35 k=2,m
35            bj(k)=bj(k)/xs
           xr=1.0d0
           xg1=bj(1)
           do 40 k=2,m
              xr=.25d0*xr*(2.0*k-3.0)**2/((k-1.0)*(2.0*k-1.0)**2)*x
40            xg1=xg1+bj(k)*xr
           xr=1.0d0
           xg2=bj(1)
           do 45 k=2,m
              xr=.25d0*xr*(2.0*k-5.0)**2/((k-1.0)*(2.0*k-3.0)**2)*x
45            xg2=xg2+bj(k)*xr
           xcs=dcos(x/2.0d0)
           xss=dsin(x/2.0d0)
           ci=el+dlog(x)-x*xss*xg1+2*xcs*xg2-2*xcs*xcs
           si=x*xcs*xg1+2*xss*xg2-dsin(x)
        else
           xr=1.0d0
           xf=1.0d0
           do 50 k=1,9
              xr=-2.0d0*xr*k*(2*k-1)/x2
50            xf=xf+xr
           xr=1.0d0/x
           xg=xr
           do 55 k=1,8
              xr=-2.0d0*xr*(2*k+1)*k/x2
55            xg=xg+xr
           ci=xf*dsin(x)/x-xg*dcos(x)/x
           si=p2-xf*dcos(x)/x-xg*dsin(x)/x
        endif
        return
        end

