        program mitairy
c
c       ===========================================================
c       purpose: this program computes the integrals of airy 
c                functions using subroutine itairy
c       input  : x   --- upper limit of the integral
c       output : apt --- integration of ai(t) from 0 and x
c                bpt --- integration of bi(t) from 0 and x
c                ant --- integration of ai(-t) from 0 and x
c                bnt --- integration of bi(-t) from 0 and x
c       example:
c
c         x      ai(t)dt       bi(t)dt       ai(-t)dt     bi(-t)dt
c        ----------------------------------------------------------
c         5    .33328759   .32147832d+03    .71788220    .15873094
c        10    .33333333   .14780980d+09    .76569840    .01504043
c        15    .33333333   .49673090d+16    .68358063    .07202621
c        20    .33333333   .47447423d+25    .71173925   -.03906173
c        25    .33333333   .78920820d+35    .70489539    .03293190
c       ===========================================================
c
        implicit double precision (a-h,o-z)
        write(*,*)'please enter x '
        read(*,*)x
        write(*,20)
        write(*,30)
        call itairy(x,apt,bpt,ant,bnt)
        write(*,10)x,apt,bpt,ant,bnt
10      format(1x,f5.1,f14.8,2x,d15.8,2f14.8)
20      format(3x,'x',8x,'ai(t)dt',7x,'bi(t)dt',9x,
     &        'ai(-t)dt',6x,'bi(-t)dt')
30      format(2x,'----------------------------------',
     &     '------------------------------')
        end


        subroutine itairy(x,apt,bpt,ant,bnt)
c
c       ======================================================
c       purpose: compute the integrals of airy fnctions with
c                respect to t from 0 and x ( x ò 0 )
c       input  : x   --- upper limit of the integral
c       output : apt --- integration of ai(t) from 0 and x
c                bpt --- integration of bi(t) from 0 and x
c                ant --- integration of ai(-t) from 0 and x
c                bnt --- integration of bi(-t) from 0 and x
c       ======================================================
c
        implicit double precision (a-h,o-z)
        dimension a(16)
        eps=1.0d-15
        pi=3.141592653589793d0
        c1=.355028053887817d0
        c2=.258819403792807d0
        sr3=1.732050807568877d0
        if (x.eq.0.0d0) then
           apt=0.0d0
           bpt=0.0d0
           ant=0.0d0
           bnt=0.0d0
        else
           if (dabs(x).le.9.25d0) then
              do 30 l=0,1
                 x=(-1)**l*x
                 fx=x
                 r=x
                 do 10 k=1,40
                    r=r*(3.0*k-2.0d0)/(3.0*k+1.0d0)*x/(3.0*k)
     &                *x/(3.0*k-1.0d0)*x
                    fx=fx+r
                    if (dabs(r).lt.dabs(fx)*eps) go to 15
10               continue
15               gx=.5d0*x*x
                 r=gx
                 do 20 k=1,40
                    r=r*(3.0*k-1.0d0)/(3.0*k+2.0d0)*x/(3.0*k)
     &                *x/(3.0*k+1.0d0)*x
                    gx=gx+r
                    if (dabs(r).lt.dabs(gx)*eps) go to 25
20               continue
25               ant=c1*fx-c2*gx
                 bnt=sr3*(c1*fx+c2*gx)
                 if (l.eq.0) then
                    apt=ant
                    bpt=bnt
                 else
                    ant=-ant
                    bnt=-bnt
                    x=-x
                 endif
30            continue
           else
              data a/.569444444444444d0,.891300154320988d0,
     &             .226624344493027d+01,.798950124766861d+01,
     &             .360688546785343d+02,.198670292131169d+03,
     &             .129223456582211d+04,.969483869669600d+04,
     &             .824184704952483d+05,.783031092490225d+06,
     &             .822210493622814d+07,.945557399360556d+08,
     &             .118195595640730d+10,.159564653040121d+11,
     &             .231369166433050d+12,.358622522796969d+13/
              q2=1.414213562373095d0
              q0=.3333333333333333d0
              q1=.6666666666666667d0
              xe=x*dsqrt(x)/1.5d0
              xp6=1.0d0/dsqrt(6.0d0*pi*xe)
              su1=1.0d0
              r=1.0d0
              xr1=1.0d0/xe
              do 35 k=1,16
                 r=-r*xr1
35               su1=su1+a(k)*r
              su2=1.0d0
              r=1.0d0
              do 40 k=1,16
                 r=r*xr1
40               su2=su2+a(k)*r
              apt=q0-dexp(-xe)*xp6*su1
              bpt=2.0d0*dexp(xe)*xp6*su2
              su3=1.0d0
              r=1.0d0
              xr2=1.0d0/(xe*xe)
              do 45 k=1,8
                 r=-r*xr2
45               su3=su3+a(2*k)*r
              su4=a(1)*xr1
              r=xr1
              do 50 k=1,7
                 r=-r*xr2
50               su4=su4+a(2*k+1)*r
              su5=su3+su4
              su6=su3-su4
              ant=q1-q2*xp6*(su5*dcos(xe)-su6*dsin(xe))
              bnt=q2*xp6*(su5*dsin(xe)+su6*dcos(xe))
           endif
        endif
        return
        end
