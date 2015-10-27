        program mairyzo
c
c       =========================================================
c       purpose: this program computes the first nt zeros of airy 
c                functions ai(x) and ai'(x), and the associated 
c                values of ai(a') and ai'(a), and the first nt 
c                zeros of airy functions bi(x) and bi'(x), and  
c                the associated values of bi(b') and bi'(b) using
c                subroutine airyzo
c       input :  nt    --- total number of zeros
c                kf    --- function code
c                          kf=1 for ai(x) and ai'(x)
c                          kf=2 for bi(x) and bi'(x)
c       output:  xa(m) --- a, the m-th zero of ai(x) or
c                          b, the m-th zero of bi(x) 
c                xb(m) --- a', the m-th zero of ai'(x) or
c                          b', the m-th zero of bi'(x)
c                xc(m) --- ai(a') or bi(b')
c                xd(m) --- ai'(a) or bi'(b)
c                          ( m --- serial number of zeros )
c       example: nt=5
c
c       m         a            ai'(a)         a'          ai(a')
c      -----------------------------------------------------------
c       1    -2.33810741     .70121082   -1.01879297    .53565666
c       2    -4.08794944    -.80311137   -3.24819758   -.41901548
c       3    -5.52055983     .86520403   -4.82009921    .38040647
c       4    -6.78670809    -.91085074   -6.16330736   -.35790794
c       5    -7.94413359     .94733571   -7.37217726    .34230124
c
c       m         b            bi'(b)         b'          bi(b')
c      -----------------------------------------------------------
c       1    -1.17371322     .60195789   -2.29443968   -.45494438
c       2    -3.27109330    -.76031014   -4.07315509    .39652284
c       3    -4.83073784     .83699101   -5.51239573   -.36796916
c       4    -6.16985213    -.88947990   -6.78129445    .34949912
c       5    -7.37676208     .92998364   -7.94017869   -.33602624
c       ==========================================================
c
        implicit double precision (a-h,o-z)
        dimension xa(50),xb(50),xc(50),xd(50)
        write(*,35)
        write(*,40)
        write(*,*)'please enter kf,nt '
        read(*,*)kf,nt
        write(*,30)kf,nt
        if (kf.eq.1) then
           write(*,*)'  m        a             ai''(a)        a''',
     &               '           ai(a'')'
        else if (kf.eq.2) then
           write(*,*)'  m        b             bi''(b)        b''',
     &               '           bi(b'')'
        endif
        write(*,*)'---------------------------------',
     &
     &            '---------------------------'
        call airyzo(nt,kf,xa,xb,xc,xd)         
        do 10 k=1,nt
10         write(*,20)k,xa(k),xd(k),xb(k),xc(k)
20      format(1x,i3,1x,3f14.8,f13.8)
30      format(1x,3hkf=,i2,',     ',3hnt=,i3)
35      format(10x,'kf=1 for ai(x) and ai''(x); kf=2 for bi(x)',
     &          ' and bi''(x)')
40      format(10x,'nt is the number of the zeros')
        end


        subroutine airyzo(nt,kf,xa,xb,xc,xd)
c
c       ========================================================
c       purpose: compute the first nt zeros of airy functions
c                ai(x) and ai'(x), a and a', and the associated
c                values of ai(a') and ai'(a); and the first nt
c                zeros of airy functions bi(x) and bi'(x), b and
c                b', and the associated values of bi(b') and
c                bi'(b)
c       input :  nt    --- total number of zeros
c                kf    --- function code
c                          kf=1 for ai(x) and ai'(x)
c                          kf=2 for bi(x) and bi'(x)
c       output:  xa(m) --- a, the m-th zero of ai(x) or
c                          b, the m-th zero of bi(x) 
c                xb(m) --- a', the m-th zero of ai'(x) or
c                          b', the m-th zero of bi'(x)
c                xc(m) --- ai(a') or bi(b')
c                xd(m) --- ai'(a) or bi'(b)
c                          ( m --- serial number of zeros )
c       routine called: airyb for computing airy functions and
c                       their derivatives
c       =======================================================
c
        implicit double precision (a-h,o-z)
        dimension xa(nt),xb(nt),xc(nt),xd(nt)
        pi=3.141592653589793d0
        do 15 i=1,nt
           if (kf.eq.1) then
              u=3.0*pi*(4.0*i-1)/8.0d0
              u1=1/(u*u)
              rt0=-(u*u)**(1.0/3.0)*((((-15.5902*u1+.929844)*u1
     &            -.138889)*u1+.10416667d0)*u1+1.0d0)
           else if (kf.eq.2) then
              if (i.eq.1) then
                 rt0=-1.17371
              else
                 u=3.0*pi*(4.0*i-3.0)/8.0
                 u1=1.0d0/(u*u)
                 rt0=-(u*u)**(1.0/3.0)*((((-15.5902*u1+.929844)*u1
     &               -.138889)*u1+.10416667)*u1+1.0)
              endif
           endif
10         x=rt0
           call airyb(x,ai,bi,ad,bd)
           if (kf.eq.1) rt=rt0-ai/ad
           if (kf.eq.2) rt=rt0-bi/bd
           if (dabs((rt-rt0)/rt).gt.1.d-9) then
              rt0=rt
              goto 10
           else
              xa(i)=rt
              if (kf.eq.1) xd(i)=ad
              if (kf.eq.2) xd(i)=bd
           endif
15      continue
        do 25 i=1,nt
           if (kf.eq.1) then
              if (i.eq.1) then
                 rt0=-1.01879
              else
                 u=3.0*pi*(4.0*i-3.0)/8.0
                 u1=1/(u*u)
                 rt0=-(u*u)**(1.0/3.0)*((((15.0168*u1-.873954)
     &            *u1+.121528)*u1-.145833d0)*u1+1.0d0)
              endif
           else if (kf.eq.2) then
              if (i.eq.1) then
                 rt0=-2.29444
              else
                 u=3.0*pi*(4.0*i-1.0)/8.0
                 u1=1.0/(u*u)
                 rt0=-(u*u)**(1.0/3.0)*((((15.0168*u1-.873954)
     &               *u1+.121528)*u1-.145833)*u1+1.0)
              endif
           endif
20         x=rt0
           call airyb(x,ai,bi,ad,bd)
           if (kf.eq.1) rt=rt0-ad/(ai*x)
           if (kf.eq.2) rt=rt0-bd/(bi*x)
           if (dabs((rt-rt0)/rt).gt.1.0d-9) then
              rt0=rt
              goto 20
           else
              xb(i)=rt
              if (kf.eq.1) xc(i)=ai
              if (kf.eq.2) xc(i)=bi
           endif
25      continue
        return
        end


        subroutine airyb(x,ai,bi,ad,bd)
c
c       =======================================================
c       purpose: compute airy functions and their derivatives
c       input:   x  --- argument of airy function
c       output:  ai --- ai(x)
c                bi --- bi(x)
c                ad --- ai'(x)
c                bd --- bi'(x)
c       =======================================================
c
        implicit double precision (a-h,o-z)
        dimension ck(41),dk(41)
        eps=1.0d-15
        pi=3.141592653589793d0
        c1=0.355028053887817d0
        c2=0.258819403792807d0
        sr3=1.732050807568877d0
        xa=dabs(x)
        xq=dsqrt(xa)
        if (x.gt.0.0d0) xm=5.0
        if (x.le.0.0d0) xm=8.0
        if (x.eq.0.0d0) then
           ai=c1
           bi=sr3*c1
           ad=-c2
           bd=sr3*c2
           return
        endif
        if (xa.le.xm) then
           fx=1.0d0
           r=1.0d0
           do 10 k=1,40
              r=r*x/(3.0d0*k)*x/(3.0d0*k-1.0d0)*x
              fx=fx+r
              if (dabs(r/fx).lt.eps) go to 15
10         continue
15         gx=x
           r=x
           do 20 k=1,40
              r=r*x/(3.0d0*k)*x/(3.0d0*k+1.0d0)*x
              gx=gx+r
              if (dabs(r/gx).lt.eps) go to 25
20         continue
25         ai=c1*fx-c2*gx
           bi=sr3*(c1*fx+c2*gx)
           df=.5d0*x*x
           r=df
           do 30 k=1,40
              r=r*x/(3.0d0*k)*x/(3.0d0*k+2.0d0)*x
              df=df+r
              if (dabs(r/df).lt.eps) go to 35
30         continue
35         dg=1.0d0
           r=1.0d0
           do 40 k=1,40
              r=r*x/(3.0d0*k)*x/(3.0d0*k-2.0d0)*x
              dg=dg+r
              if (dabs(r/dg).lt.eps) go to 45
40         continue
45         ad=c1*df-c2*dg
           bd=sr3*(c1*df+c2*dg)
        else
           xe=xa*xq/1.5d0
           xr1=1.0d0/xe
           xar=1.0d0/xq
           xf=dsqrt(xar)
           rp=.5641895835477563d0
           r=1.0d0
           do 50 k=1,40
              r=r*(6.0d0*k-1.0d0)/216.0d0*(6.0d0*k-3.0d0)
     &          /k*(6.0d0*k-5.0d0)/(2.0d0*k-1.0d0)
              ck(k)=r
50            dk(k)=-(6.0d0*k+1.0d0)/(6.0d0*k-1.0d0)*ck(k)
           km=int(24.5-xa)
           if (xa.lt.6.0) km=14
           if (xa.gt.15.0) km=10
           if (x.gt.0.0d0) then
              sai=1.0d0
              sad=1.0d0
              r=1.0d0
              do 55 k=1,km
                 r=-r*xr1
                 sai=sai+ck(k)*r
55               sad=sad+dk(k)*r
              sbi=1.0d0
              sbd=1.0d0
              r=1.0d0
              do 60 k=1,km
                 r=r*xr1
                 sbi=sbi+ck(k)*r
60               sbd=sbd+dk(k)*r
              xp1=dexp(-xe)
              ai=.5d0*rp*xf*xp1*sai
              bi=rp*xf/xp1*sbi
              ad=-.5d0*rp/xf*xp1*sad
              bd=rp/xf/xp1*sbd
           else
              xcs=dcos(xe+pi/4.0d0)
              xss=dsin(xe+pi/4.0d0)
              ssa=1.0d0
              sda=1.0d0
              r=1.0d0
              xr2=1.0d0/(xe*xe)
              do 65 k=1,km
                 r=-r*xr2
                 ssa=ssa+ck(2*k)*r
65               sda=sda+dk(2*k)*r
              ssb=ck(1)*xr1
              sdb=dk(1)*xr1
              r=xr1
              do 70 k=1,km
                 r=-r*xr2
                 ssb=ssb+ck(2*k+1)*r
70               sdb=sdb+dk(2*k+1)*r
              ai=rp*xf*(xss*ssa-xcs*ssb)
              bi=rp*xf*(xcs*ssa+xss*ssb)
              ad=-rp/xf*(xcs*sda+xss*sdb)
              bd=rp/xf*(xss*sda-xcs*sdb)
           endif
        endif
        return
        end
