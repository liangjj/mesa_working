        program mairya
c
c       ============================================================
c       purpose: this program computes airy functions and their 
c                derivatives using subroutine airya
c       input:   x  --- argument of airy function
c       output:  ai --- ai(x)
c                bi --- bi(x)
c                ad --- ai'(x)
c                bd --- bi'(x)
c       example:
c
c   x       ai(x)          bi(x)          ai'(x)         bi'(x)
c  ----------------------------------------------------------------
c   0   .35502805d+00  .61492663d+00 -.25881940d+00  .44828836d+00
c  10   .11047533d-09  .45564115d+09 -.35206337d-09  .14292361d+10
c  20   .16916729d-26  .21037650d+26 -.75863916d-26  .93818393d+26
c  30   .32082176d-48  .90572885d+47 -.17598766d-47  .49533045d+48
c
c   x       ai(-x)         bi(-x)         ai'(-x)        bi'(-x)
c  ----------------------------------------------------------------
c   0       .35502805      .61492663     -.25881940      .44828836
c  10       .04024124     -.31467983      .99626504      .11941411
c  20      -.17640613     -.20013931      .89286286     -.79142903
c  30      -.08796819     -.22444694     1.22862060     -.48369473
c       ============================================================
c
        implicit double precision (a-h,o-z)
        write(*,*)'please enter x '
        read(*,*)x
        call airya(x,ai,bi,ad,bd)
        write(*,30)
        write(*,40)
        write(*,10)x,ai,bi,ad,bd
        write(*,*)
        call airya(-x,ai,bi,ad,bd)
        write(*,50)
        write(*,40)
        write(*,20)x,ai,bi,ad,bd
10      format(1x,f5.1,4d16.8)
20      format(1x,f5.1,4d16.8)
30      format(4x,'x',8x,'ai(x)',11x,'bi(x)',11x,'ai''(x)',
     &         10x,'bi''(x)')
40      format(2x,'----------------------------------',
     &        '-----------------------------------')
50      format(4x,'x',8x,'ai(-x)',10x,'bi(-x)',10x,
     &        'ai''(-x)',9x,'bi''(-x)')
        end


        subroutine airya(x,ai,bi,ad,bd)
c
c       ======================================================
c       purpose: compute airy functions and their derivatives
c       input:   x  --- argument of airy function
c       output:  ai --- ai(x)
c                bi --- bi(x)
c                ad --- ai'(x)
c                bd --- bi'(x)
c       routine called:
c                ajyik for computing jv(x), yv(x), iv(x) and
c                kv(x) with v=1/3 and 2/3
c       ======================================================
c
        implicit double precision (a-h,o-z)
        xa=dabs(x)
        pir=0.318309886183891d0
        c1=0.355028053887817d0
        c2=0.258819403792807d0
        sr3=1.732050807568877d0
        z=xa**1.5/1.5d0
        xq=dsqrt(xa)
        call ajyik(z,vj1,vj2,vy1,vy2,vi1,vi2,vk1,vk2)
        if (x.eq.0.0d0) then
           ai=c1
           bi=sr3*c1
           ad=-c2
           bd=sr3*c2
        else if (x.gt.0.0d0) then
           ai=pir*xq/sr3*vk1
           bi=xq*(pir*vk1+2.0d0/sr3*vi1)
           ad=-xa/sr3*pir*vk2
           bd=xa*(pir*vk2+2.0d0/sr3*vi2)
        else
           ai=0.5d0*xq*(vj1-vy1/sr3)
           bi=-0.5d0*xq*(vj1/sr3+vy1)
           ad=0.5d0*xa*(vj2+vy2/sr3)
           bd=0.5d0*xa*(vj2/sr3-vy2)
        endif
        return
        end


        subroutine ajyik(x,vj1,vj2,vy1,vy2,vi1,vi2,vk1,vk2)
c
c       =======================================================
c       purpose: compute bessel functions jv(x) and yv(x),
c                and modified bessel functions iv(x) and
c                kv(x), and their derivatives with v=1/3,2/3
c       input :  x --- argument of jv(x),yv(x),iv(x) and
c                      kv(x) ( x ò 0 )
c       output:  vj1 --- j1/3(x)
c                vj2 --- j2/3(x)
c                vy1 --- y1/3(x)
c                vy2 --- y2/3(x)
c                vi1 --- i1/3(x)
c                vi2 --- i2/3(x)
c                vk1 --- k1/3(x)
c                vk2 --- k2/3(x)
c       =======================================================
c
        implicit double precision (a-h,o-z)
        if (x.eq.0.0d0) then
           vj1=0.0d0
           vj2=0.0d0
           vy1=-1.0d+300
           vy2=1.0d+300
           vi1=0.0d0
           vi2=0.0d0
           vk1=-1.0d+300
           vk2=-1.0d+300
           return
        endif
        pi=3.141592653589793d0
        rp2=.63661977236758d0
        gp1=.892979511569249d0
        gp2=.902745292950934d0
        gn1=1.3541179394264d0
        gn2=2.678938534707747d0
        vv0=0.444444444444444d0
        uu0=1.1547005383793d0
        x2=x*x
        k0=12
        if (x.ge.35.0) k0=10
        if (x.ge.50.0) k0=8
        if (x.le.12.0) then
           do 25 l=1,2
              vl=l/3.0d0
              vjl=1.0d0
              r=1.0d0
              do 15 k=1,40
                 r=-0.25d0*r*x2/(k*(k+vl))
                 vjl=vjl+r
                 if (dabs(r).lt.1.0d-15) go to 20
15            continue
20            a0=(0.5d0*x)**vl
              if (l.eq.1) vj1=a0/gp1*vjl
              if (l.eq.2) vj2=a0/gp2*vjl
25         continue
        else
           do 40 l=1,2
              vv=vv0*l*l
              px=1.0d0
              rp=1.0d0
              do 30 k=1,k0
                 rp=-0.78125d-2*rp*(vv-(4.0*k-3.0)**2.0)*(vv-
     &              (4.0*k-1.0)**2.0)/(k*(2.0*k-1.0)*x2)
30               px=px+rp
              qx=1.0d0
              rq=1.0d0
              do 35 k=1,k0
                 rq=-0.78125d-2*rq*(vv-(4.0*k-1.0)**2.0)*(vv-
     &              (4.0*k+1.0)**2.0)/(k*(2.0*k+1.0)*x2)
35               qx=qx+rq
              qx=0.125d0*(vv-1.0)*qx/x
              xk=x-(0.5d0*l/3.0d0+0.25d0)*pi
              a0=dsqrt(rp2/x)
              ck=dcos(xk)
              sk=dsin(xk)
              if (l.eq.1) then
                 vj1=a0*(px*ck-qx*sk)
                 vy1=a0*(px*sk+qx*ck)
              else if (l.eq.2) then
                 vj2=a0*(px*ck-qx*sk)
                 vy2=a0*(px*sk+qx*ck)
              endif
40         continue
        endif
        if (x.le.12.0d0) then
           do 55 l=1,2
              vl=l/3.0d0
              vjl=1.0d0
              r=1.0d0
              do 45 k=1,40
                 r=-0.25d0*r*x2/(k*(k-vl))
                 vjl=vjl+r
                 if (dabs(r).lt.1.0d-15) go to 50
45            continue
50            b0=(2.0d0/x)**vl
              if (l.eq.1) uj1=b0*vjl/gn1
              if (l.eq.2) uj2=b0*vjl/gn2
55         continue
           pv1=pi/3.0d0
           pv2=pi/1.5d0
           vy1=uu0*(vj1*dcos(pv1)-uj1)
           vy2=uu0*(vj2*dcos(pv2)-uj2)
        endif
        if (x.le.18.0) then
           do 70 l=1,2
              vl=l/3.0d0
              vil=1.0d0
              r=1.0d0
              do 60 k=1,40
                 r=0.25d0*r*x2/(k*(k+vl))
                 vil=vil+r
                 if (dabs(r).lt.1.0d-15) go to 65
60            continue
65            a0=(0.5d0*x)**vl
              if (l.eq.1) vi1=a0/gp1*vil
              if (l.eq.2) vi2=a0/gp2*vil
70         continue
        else
           c0=dexp(x)/dsqrt(2.0d0*pi*x)
           do 80 l=1,2
              vv=vv0*l*l
              vsl=1.0d0
              r=1.0d0
              do 75 k=1,k0
                 r=-0.125d0*r*(vv-(2.0d0*k-1.0d0)**2.0)/(k*x)
75               vsl=vsl+r
              if (l.eq.1) vi1=c0*vsl
              if (l.eq.2) vi2=c0*vsl
80         continue
        endif
        if (x.le.9.0d0) then
           do 95 l=1,2
              vl=l/3.0d0
               if (l.eq.1) gn=gn1
               if (l.eq.2) gn=gn2
               a0=(2.0d0/x)**vl/gn
               sum=1.0d0
               r=1.0d0
               do 85 k=1,60
                  r=0.25d0*r*x2/(k*(k-vl))
                  sum=sum+r
                  if (dabs(r).lt.1.0d-15) go to 90
85             continue
90            if (l.eq.1) vk1=0.5d0*uu0*pi*(sum*a0-vi1)
              if (l.eq.2) vk2=0.5d0*uu0*pi*(sum*a0-vi2)
95         continue
        else
           c0=dexp(-x)*dsqrt(0.5d0*pi/x)
           do 105 l=1,2
              vv=vv0*l*l
              sum=1.0d0
              r=1.0d0
              do 100 k=1,k0
                 r=0.125d0*r*(vv-(2.0*k-1.0)**2.0)/(k*x)
100              sum=sum+r
              if (l.eq.1) vk1=c0*sum
              if (l.eq.2) vk2=c0*sum
105        continue
        endif
        return
        end


