        program mklvnb
c
c       ========================================================
c       purpose: this program computes kelvin functions ber x, 
c                bei x, ker x and kei x, and their derivatives 
c                using subroutine klvnb
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
c     x       ber x         bei x         ker x         kei x
c   -------------------------------------------------------------
c     0    .100000d+01   .000000d+00    ì           -.785398d+00
c     5   -.623008d+01   .116034d+00  -.115117d-01   .111876d-01
c    10    .138840d+03   .563705d+02   .129466d-03  -.307525d-03
c    15   -.296725d+04  -.295271d+04  -.151433d-07   .796289d-05
c    20    .474894d+05   .114775d+06  -.771523d-07  -.185894d-06
c
c     x       ber'x         bei'x         ker'x         kei'x
c   -------------------------------------------------------------
c     0    .000000d+00   .000000d+00  - ì            .000000d+00
c     5   -.384534d+01  -.435414d+01   .171934d-01  -.819979d-03
c    10    .511952d+02   .135309d+03  -.315597d-03   .140914d-03
c    15    .910555d+02  -.408776d+04   .564468d-05  -.588222d-05
c    20   -.488032d+05   .111855d+06  -.750186d-07   .190624d-06
c       ========================================================
c
        implicit double precision (a-h,o-z)
        write(*,*)'please enter x '
        read(*,*)x
        write(*,*)'   x        ber x           bei x',
     &            '           ker x           kei x'
        write(*,*)'--------------------------------',
     &            '--------------------------------------'
        call klvnb(x,ber,bei,ger,gei,der,dei,her,hei)
        write(*,10)x,ber,bei,ger,gei
        write(*,*)
        write(*,*)'   x        ber''x           bei''x',
     &            '           ker''x           kei''x'
        write(*,*)'--------------------------------',
     &            '--------------------------------------'
        write(*,10)x,der,dei,her,hei
10      format(1x,f5.1,4d16.6)
        end


        subroutine klvnb(x,ber,bei,ger,gei,der,dei,her,hei)
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
        if (x.eq.0.0d0) then
           ber=1.0d0
           bei=0.0d0
           ger=1.0d+300
           gei=-.25d0*pi
           der=0.0d0
           dei=0.0d0
           her=-1.0d+300
           hei=0.0d0
        else if (x.lt.8.0d0) then
           t=x/8.0d0
           t2=t*t
           u=t2*t2
           ber=((((((-.901d-5*u+.122552d-2)*u-.08349609d0)*u
     &         +2.64191397d0)*u-32.36345652d0)*u
     &         +113.77777774d0)*u-64.0d0)*u+1.0d0
           bei=t*t*((((((.11346d-3*u-.01103667d0)*u
     &         +.52185615d0)*u-10.56765779d0)*u
     &         +72.81777742d0)*u-113.77777774d0)*u+16.0d0)
           ger=((((((-.2458d-4*u+.309699d-2)*u-.19636347d0)
     &         *u+5.65539121d0)*u-60.60977451d0)*u+
     &         171.36272133d0)*u-59.05819744d0)*u-.57721566d0
           ger=ger-dlog(.5d0*x)*ber+.25d0*pi*bei
           gei=t2*((((((.29532d-3*u-.02695875d0)*u
     &         +1.17509064d0)*u-21.30060904d0)*u
     &         +124.2356965d0)*u-142.91827687d0)*u
     &         +6.76454936d0)
           gei=gei-dlog(.5d0*x)*bei-.25d0*pi*ber
           der=x*t2*((((((-.394d-5*u+.45957d-3)*u
     &         -.02609253d0)*u+.66047849d0)*u-6.0681481d0)*u
     &         +14.22222222d0)*u-4.0d0)
           dei=x*((((((.4609d-4*u-.379386d-2)*u+.14677204d0)
     &         *u-2.31167514d0)*u+11.37777772d0)*u
     &         -10.66666666d0)*u+.5d0)
           her=x*t2*((((((-.1075d-4*u+.116137d-2)*u
     &         -.06136358d0)*u+1.4138478d0)*u-11.36433272d0)
     &         *u+21.42034017d0)*u-3.69113734d0)
           her=her-dlog(.5d0*x)*der-ber/x+.25d0*pi*dei
           hei=x*((((((.11997d-3*u-.926707d-2)*u
     &         +.33049424d0)*u-4.65950823d0)*u+19.41182758d0)
     &         *u-13.39858846d0)*u+.21139217d0)
           hei=hei-dlog(.5d0*x)*dei-bei/x-.25d0*pi*der
        else
           t=8.0d0/x
           do 10 l=1,2
              v=(-1)**l*t
              tpr=((((.6d-6*v-.34d-5)*v-.252d-4)*v-.906d-4)
     &            *v*v+.0110486d0)*v
              tpi=((((.19d-5*v+.51d-5)*v*v-.901d-4)*v
     &            -.9765d-3)*v-.0110485d0)*v-.3926991d0
              if (l.eq.1) then
                 tnr=tpr
                 tni=tpi
              endif
10         continue
           yd=x/dsqrt(2.0d0)
           ye1=dexp(yd+tpr)
           ye2=dexp(-yd+tnr)
           yc1=1.0d0/dsqrt(2.0d0*pi*x)
           yc2=dsqrt(pi/(2.0d0*x))
           csp=dcos(yd+tpi)
           ssp=dsin(yd+tpi)
           csn=dcos(-yd+tni)
           ssn=dsin(-yd+tni)
           ger=yc2*ye2*csn
           gei=yc2*ye2*ssn
           fxr=yc1*ye1*csp
           fxi=yc1*ye1*ssp
           ber=fxr-gei/pi
           bei=fxi+ger/pi
           do 15 l=1,2
              v=(-1)**l*t
              ppr=(((((.16d-5*v+.117d-4)*v+.346d-4)*v+.5d-6)
     &            *v-.13813d-2)*v-.0625001d0)*v+.7071068d0
              ppi=(((((-.32d-5*v-.24d-5)*v+.338d-4)*v+
     &           .2452d-3)*v+.13811d-2)*v-.1d-6)*v+.7071068d0
              if (l.eq.1) then
                 pnr=ppr
                 pni=ppi
              endif
15         continue
           her=gei*pni-ger*pnr
           hei=-(gei*pnr+ger*pni)
           der=fxr*ppr-fxi*ppi-hei/pi
           dei=fxi*ppr+fxr*ppi+her/pi
        endif
        return
        end
