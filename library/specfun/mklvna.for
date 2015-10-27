        program mklvna
c
c       =======================================================
c       purpose: this program computes kelvin functions ber x, 
c                bei x, ker x and kei x, and their derivatives 
c                using subroutine klvna
c       input :  x   --- argument of kelvin functions
c       output:  ber --- ber x
c                bei --- bei x
c                ger --- ker x
c                gei --- kei x
c                der --- ber'x
c                dei --- bei'x
c                her --- ker'x
c                hei --- kei'x
c       example:
c
c      x       ber x          bei x          ker x          kei x
c    -----------------------------------------------------------------
c      0    .1000000d+01    0              ì            -.7853982d+00
c      5   -.6230082d+01   .1160344d+00  -.1151173d-01   .1118759d-01
c     10    .1388405d+03   .5637046d+02   .1294663d-03  -.3075246d-03
c     15   -.2967255d+04  -.2952708d+04  -.1514347d-07   .7962894d-05
c     20    .4748937d+05   .1147752d+06  -.7715233d-07  -.1858942d-06
c
c      x       ber'x          bei'x          ker'x          kei'x
c    -----------------------------------------------------------------
c      0     0              0            - ì              0
c      5   -.3845339d+01  -.4354141d+01   .1719340d-01  -.8199865d-03
c     10    .5119526d+02   .1353093d+03  -.3155969d-03   .1409138d-03
c     15    .9105533d+02  -.4087755d+04   .5644678d-05  -.5882223d-05
c     20   -.4880320d+05   .1118550d+06  -.7501859d-07   .1906243d-06
c       =======================================================
c
        implicit double precision (a-h,o-z)
        write(*,*)'please enter x '
        read(*,*)x
        call klvna(x,ber,bei,ger,gei,der,dei,her,hei)
        write(*,*)'   x        ber x           bei x',
     &            '           ker x           kei x'
        write(*,*)'--------------------------------',
     &            '--------------------------------------'
        write(*,10)x,ber,bei,ger,gei
        write(*,*)
        write(*,*)'   x        ber''x           bei''x',
     &            '           ker''x           kei''x'
        write(*,*)'--------------------------------',
     &            '--------------------------------------'
        write(*,10)x,der,dei,her,hei
10      format(1x,f5.1,4d16.8)
        end


        subroutine klvna(x,ber,bei,ger,gei,der,dei,her,hei)
c
c       ======================================================
c       purpose: compute kelvin functions ber x, bei x, ker x
c                and kei x, and their derivatives  ( x > 0 )
c       input :  x   --- argument of kelvin functions
c       output:  ber --- ber x
c                bei --- bei x
c                ger --- ker x
c                gei --- kei x
c                der --- ber'x
c                dei --- bei'x
c                her --- ker'x
c                hei --- kei'x
c       ================================================
c
        implicit double precision (a-h,o-z)
        pi=3.141592653589793d0
        el=.5772156649015329d0
        eps=1.0d-15
        if (x.eq.0.0d0) then
           ber=1.0d0
           bei=0.0d0
           ger=1.0d+300
           gei=-0.25d0*pi
           der=0.0d0
           dei=0.0d0
           her=-1.0d+300
           hei=0.0d0
           return
        endif
        x2=0.25d0*x*x
        x4=x2*x2
        if (dabs(x).lt.10.0d0) then
           ber=1.0d0
           r=1.0d0
           do 10 m=1,60
              r=-0.25d0*r/(m*m)/(2.0d0*m-1.0d0)**2*x4
              ber=ber+r
              if (dabs(r).lt.dabs(ber)*eps) go to 15
10         continue
15         bei=x2
           r=x2
           do 20 m=1,60
              r=-0.25d0*r/(m*m)/(2.0d0*m+1.0d0)**2*x4
              bei=bei+r
              if (dabs(r).lt.dabs(bei)*eps) go to 25
20         continue
25         ger=-(dlog(x/2.0d0)+el)*ber+0.25d0*pi*bei
           r=1.0d0
           gs=0.0d0
           do 30 m=1,60
              r=-0.25d0*r/(m*m)/(2.0d0*m-1.0d0)**2*x4
              gs=gs+1.0d0/(2.0d0*m-1.0d0)+1.0d0/(2.0d0*m)
              ger=ger+r*gs
              if (dabs(r*gs).lt.dabs(ger)*eps) go to 35
30         continue
35         gei=x2-(dlog(x/2.0d0)+el)*bei-0.25d0*pi*ber
           r=x2
           gs=1.0d0
           do 40 m=1,60
              r=-0.25d0*r/(m*m)/(2.0d0*m+1.0d0)**2*x4
              gs=gs+1.0d0/(2.0d0*m)+1.0d0/(2.0d0*m+1.0d0)
              gei=gei+r*gs
              if (dabs(r*gs).lt.dabs(gei)*eps) go to 45
40         continue
45         der=-0.25d0*x*x2
           r=der
           do 50 m=1,60
              r=-0.25d0*r/m/(m+1.0d0)/(2.0d0*m+1.0d0)**2*x4
              der=der+r
              if (dabs(r).lt.dabs(der)*eps) go to 55
50         continue
55         dei=0.5d0*x
           r=dei
           do 60 m=1,60
              r=-0.25d0*r/(m*m)/(2.d0*m-1.d0)/(2.d0*m+1.d0)*x4
              dei=dei+r
              if (dabs(r).lt.dabs(dei)*eps) go to 65
60            continue
65         r=-0.25d0*x*x2
           gs=1.5d0
           her=1.5d0*r-ber/x-(dlog(x/2.d0)+el)*der+0.25*pi*dei
           do 70 m=1,60
              r=-0.25d0*r/m/(m+1.0d0)/(2.0d0*m+1.0d0)**2*x4
              gs=gs+1.0d0/(2*m+1.0d0)+1.0d0/(2*m+2.0d0)
              her=her+r*gs
              if (dabs(r*gs).lt.dabs(her)*eps) go to 75
70         continue
75         r=0.5d0*x
           gs=1.0d0
           hei=0.5d0*x-bei/x-(dlog(x/2.d0)+el)*dei-0.25*pi*der
           do 80 m=1,60
              r=-0.25d0*r/(m*m)/(2*m-1.0d0)/(2*m+1.0d0)*x4
              gs=gs+1.0d0/(2.0d0*m)+1.0d0/(2*m+1.0d0)
              hei=hei+r*gs
              if (dabs(r*gs).lt.dabs(hei)*eps) return
80         continue
        else
           pp0=1.0d0
           pn0=1.0d0
           qp0=0.0d0
           qn0=0.0d0
           r0=1.0d0
           km=18
           if (dabs(x).ge.40.0) km=10
           fac=1.0d0
           do 85 k=1,km
              fac=-fac
              xt=0.25d0*k*pi-int(0.125d0*k)*2.0d0*pi
              cs=cos(xt)
              ss=sin(xt)
              r0=0.125d0*r0*(2.0d0*k-1.0d0)**2/k/x
              rc=r0*cs
              rs=r0*ss
              pp0=pp0+rc
              pn0=pn0+fac*rc
              qp0=qp0+rs
85            qn0=qn0+fac*rs
           xd=x/dsqrt(2.0d0)
           xe1=dexp(xd)
           xe2=dexp(-xd)
           xc1=1.d0/dsqrt(2.0d0*pi*x)
           xc2=dsqrt(.5d0*pi/x)
           cp0=dcos(xd+0.125d0*pi)
           cn0=dcos(xd-0.125d0*pi)
           sp0=dsin(xd+0.125d0*pi)
           sn0=dsin(xd-0.125d0*pi)
           ger=xc2*xe2*(pn0*cp0-qn0*sp0)
           gei=xc2*xe2*(-pn0*sp0-qn0*cp0)
           ber=xc1*xe1*(pp0*cn0+qp0*sn0)-gei/pi
           bei=xc1*xe1*(pp0*sn0-qp0*cn0)+ger/pi
           pp1=1.0d0
           pn1=1.0d0
           qp1=0.0d0
           qn1=0.0d0
           r1=1.0d0
           fac=1.0d0
           do 90 k=1,km
              fac=-fac
              xt=0.25d0*k*pi-int(0.125d0*k)*2.0d0*pi
              cs=dcos(xt)
              ss=dsin(xt)
              r1=0.125d0*r1*(4.d0-(2.0d0*k-1.0d0)**2)/k/x
              rc=r1*cs
              rs=r1*ss
              pp1=pp1+fac*rc
              pn1=pn1+rc
              qp1=qp1+fac*rs
              qn1=qn1+rs
90         continue
           her=xc2*xe2*(-pn1*cn0+qn1*sn0)
           hei=xc2*xe2*(pn1*sn0+qn1*cn0)
           der=xc1*xe1*(pp1*cp0+qp1*sp0)-hei/pi
           dei=xc1*xe1*(pp1*sp0-qp1*cp0)+her/pi
        endif
        return
        end
