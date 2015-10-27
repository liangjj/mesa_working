        program mklvnzo
c
c       ==============================================================
c       purpose: this program computes the first nt zeros of kelvin 
c                functions and their derivatives using subroutine
c                klvnzo
c       input :  nt --- total number of zeros
c       example: nt = 5
c
c           zeros of kelvin functions ber x, bei x, ker x and kei x
c
c        m       ber x          bei x          ker x          kei x
c       ---------------------------------------------------------------
c        1     2.84891782     5.02622395     1.71854296     3.91466761
c        2     7.23882945     9.45540630     6.12727913     8.34422506
c        3    11.67396355    13.89348785    10.56294271    12.78255715
c        4    16.11356383    18.33398346    15.00268812    17.22314372
c        5    20.55463158    22.77543929    19.44381663    21.66464214
c
c          zeros of kelvin functions ber'x, bei'x, ker'x and kei'x
c
c        m       ber'x          bei'x          ker'x          kei'x
c       ---------------------------------------------------------------
c        1     6.03871081     3.77267330     2.66583979     4.93181194
c        2    10.51364251     8.28098785     7.17212212     9.40405458
c        3    14.96844542    12.74214752    11.63218639    13.85826916
c        4    19.41757493    17.19343175    16.08312025    18.30717294
c        5    23.86430432    21.64114394    20.53067845    22.75379258
c       ==============================================================
c
        implicit double precision (a-h,o-z)
        dimension r1(50),r2(50),r3(50),r4(50),r5(50),r6(50),
     &            r7(50),r8(50)
        write(*,35)
        write(*,*)'please enter nt '
        read(*,*)nt
        write(*,25)
        write(*,*)
        write(*,*)'   m       ber x          bei x          ker x',
     &            '          kei x'
        write(*,*)' ------------------------------------------------',
     &            '---------------'
        call klvnzo(nt,1,r1)
        call klvnzo(nt,2,r2)
        call klvnzo(nt,3,r3)
        call klvnzo(nt,4,r4)
        do 10 l=1,nt
10         write(*,20)l,r1(l),r2(l),r3(l),r4(l)
        call klvnzo(nt,5,r5)
        call klvnzo(nt,6,r6)
        call klvnzo(nt,7,r7)
        call klvnzo(nt,8,r8)
        write(*,*)
        write(*,30)
        write(*,*)
        write(*,*)'   m       ber''x          bei''x          ker''x',
     &            '          kei''x'
        write(*,*)' ------------------------------------------------',
     &            '---------------'
        do 15 l=1,nt
15         write(*,20)l,r5(l),r6(l),r7(l),r8(l)
        close (01)
20      format(1x,i3,1x,f14.8,1x,f14.8,1x,f14.8,1x,f14.8)
25      format(4x,'zeros of kelvin functions ber x, bei x,'
     &        ,' ker x and kei x')
30      format(4x,'zeros of kelvin functions ber''x, bei''x,'
     &        ,' ker''x and kei''x')
35      format(1x,'nt is the number of the zeros')
        end


        subroutine klvnzo(nt,kd,zo)
c
c       ====================================================
c       purpose: compute the zeros of kelvin functions
c       input :  nt  --- total number of zeros
c                kd  --- function code
c                kd=1 to 8 for ber x, bei x, ker x, kei x,
c                          ber'x, bei'x, ker'x and kei'x,
c                          respectively.
c       output:  zo(m) --- the m-th zero of kelvin function
c                          for code kd
c       routine called:
c                klvna for computing kelvin functions and
c                their derivatives
c       ====================================================
c
        implicit double precision (a-h,o-z)
        dimension zo(nt),rt0(8)
        rt0(1)=2.84891
        rt0(2)=5.02622
        rt0(3)=1.71854
        rt0(4)=3.91467
        rt0(5)=6.03871
        rt0(6)=3.77268
        rt0(7)=2.66584
        rt0(8)=4.93181
        rt=rt0(kd)
        do 15 m=1,nt
10         call klvna(rt,ber,bei,ger,gei,der,dei,her,hei)
           if (kd.eq.1) then
              rt=rt-ber/der
           else if (kd.eq.2) then
              rt=rt-bei/dei
           else if (kd.eq.3) then
              rt=rt-ger/her
           else if (kd.eq.4) then
              rt=rt-gei/hei
           else if (kd.eq.5) then
              ddr=-bei-der/rt
              rt=rt-der/ddr
           else if (kd.eq.6) then
              ddi=ber-dei/rt
              rt=rt-dei/ddi
           else if (kd.eq.7) then
              gdr=-gei-her/rt
              rt=rt-her/gdr
           else
              gdi=ger-hei/rt
              rt=rt-hei/gdi
           endif
           if (dabs(rt-rt0(kd)).gt.5.0d-10) then
              rt0(kd)=rt
              go to 10
           endif
           zo(m)=rt
15         rt=rt+4.44d0
        return
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
           gei=-.25d0*pi
           der=0.0d0
           dei=0.0d0
           her=-1.0d+300
           hei=0.0d0
           return
        endif
        x2=.25d0*x*x
        x4=x2*x2
        if (dabs(x).lt.10.0d0) then
           ber=1.0d0
           r=1.0d0
           do 10 m=1,60
              r=-.25d0*r/(m*m)/(2.0d0*m-1.0d0)**2*x4
              ber=ber+r
              if (dabs(r/ber).lt.eps) go to 15
10         continue
15         bei=x2
           r=x2
           do 20 m=1,60
              r=-.25d0*r/(m*m)/(2.0d0*m+1.0d0)**2*x4
              bei=bei+r
              if (dabs(r/bei).lt.eps) go to 25
20         continue
25         ger=-(dlog(x/2.0d0)+el)*ber+.25d0*pi*bei
           r=1.0d0
           gs=0.0d0
           do 30 m=1,60
              r=-.25d0*r/(m*m)/(2.0d0*m-1.0d0)**2*x4
              gs=gs+1.0d0/(2.0d0*m-1.0d0)+1.0d0/(2.0d0*m)
              ger=ger+r*gs
              if (dabs(r*gs/ger).lt.eps) go to 35
30         continue
35         gei=x2-(dlog(x/2.0d0)+el)*bei-.25d0*pi*ber
           r=x2
           gs=1.0d0
           do 40 m=1,60
              r=-.25d0*r/(m*m)/(2.0d0*m+1.0d0)**2*x4
              gs=gs+1.0d0/(2.0d0*m)+1.0d0/(2.0d0*m+1.0d0)
              gei=gei+r*gs
              if (dabs(r*gs/gei).lt.eps) go to 45
40         continue
45         der=-.25d0*x*x2
           r=der
           do 50 m=1,60
              r=-.25d0*r/m/(m+1.0d0)/(2.0d0*m+1.0d0)**2*x4
              der=der+r
              if (dabs(r/der).lt.eps) go to 55
50         continue
55         dei=.5d0*x
           r=dei
           do 60 m=1,60
              r=-.25d0*r/(m*m)/(2.d0*m-1.d0)/(2.d0*m+1.d0)*x4
              dei=dei+r
              if (dabs(r/dei).lt.eps) go to 65
60            continue
65         r=-.25d0*x*x2
           gs=1.5d0
           her=1.5d0*r-ber/x-(dlog(x/2.d0)+el)*der+.25*pi*dei
           do 70 m=1,60
              r=-.25d0*r/m/(m+1.0d0)/(2.0d0*m+1.0d0)**2*x4
              gs=gs+1.0d0/(2*m+1.0d0)+1.0d0/(2*m+2.0d0)
              her=her+r*gs
              if (dabs(r*gs/her).lt.eps) go to 75
70         continue
75         r=.5d0*x
           gs=1.0d0
           hei=.5d0*x-bei/x-(dlog(x/2.d0)+el)*dei-.25*pi*der
           do 80 m=1,60
              r=-.25d0*r/(m*m)/(2*m-1.0d0)/(2*m+1.0d0)*x4
              gs=gs+1.0d0/(2.0d0*m)+1.0d0/(2*m+1.0d0)
              hei=hei+r*gs
              if (dabs(r*gs/hei).lt.eps) return
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
              xt=.25d0*k*pi-int(.125d0*k)*2.0d0*pi
              cs=cos(xt)
              ss=sin(xt)
              r0=.125d0*r0*(2.0d0*k-1.0d0)**2/k/x
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
           cp0=dcos(xd+.125d0*pi)
           cn0=dcos(xd-.125d0*pi)
           sp0=dsin(xd+.125d0*pi)
           sn0=dsin(xd-.125d0*pi)
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
              xt=.25d0*k*pi-int(.125d0*k)*2.0d0*pi
              cs=dcos(xt)
              ss=dsin(xt)
              r1=.125d0*r1*(4.d0-(2.0d0*k-1.0d0)**2)/k/x
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
