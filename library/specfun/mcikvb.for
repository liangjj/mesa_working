        program mcikvb
c
c       =============================================================
c       purpose: this program computes the modified bessel functions 
c                iv(z), kv(z) and their derivatives for an arbitrary
c                order and complex argument using subroutine cikvb 
c       input :  z --- complex argument z
c                v --- real order of iv(z) and kv(z)
c                      ( v = n+v0,  0 ף n ף 250, 0 ף v0 < 1 )
c       output:  cbi(n) --- in+v0(z)
c                cdi(n) --- in+v0'(z)
c                cbk(n) --- kn+v0(z)
c                cdk(n) --- kn+v0'(z)
c       example: v= n+v0,   v0 = .25,   z = 4.0+ i 2.0
c
c    n     re[iv(z)]      im[iv(z)]     re[iv'(z)]     im[iv'(z)]
c  -----------------------------------------------------------------
c    0  -.19336550d+01  .10328998d+02 -.23119621d+01  .91612230d+01
c    1  -.24735044d+01  .85964317d+01 -.23898329d+01  .78707023d+01
c    2  -.28460107d+01  .54124063d+01 -.24105909d+01  .55204965d+01
c    3  -.23476775d+01  .24445612d+01 -.21145027d+01  .30604463d+01
c    4  -.13829947d+01  .70848630d+00 -.14732387d+01  .12545751d+01
c    5  -.59879982d+00  .64588999d-01 -.78816416d+00  .32629794d+00
c
c    n     re[kv(z)]      im[kv(z)]     re[kv'(z)]     im[kv'(z)]
c   ----------------------------------------------------------------
c    0  -.64820386d-02 -.84715754d-02  .75118612d-02  .89920077d-02
c    1  -.80477525d-02 -.92535355d-02  .96506687d-02  .97789903d-02
c    2  -.12819299d-01 -.11086405d-01  .16310878d-01  .11358076d-01
c    3  -.24574004d-01 -.13462616d-01  .33167751d-01  .11850554d-01
c    4  -.53516204d-01 -.12614703d-01  .75424026d-01  .14407268d-02
c    5  -.12627405d+00  .10581162d-01  .18054884d+00 -.64789392d-01
c       =============================================================
c
        implicit double precision (a,b,d-h,o-y)
        implicit complex*16 (c,z)
        dimension cbi(0:250),cdi(0:250),cbk(0:250),cdk(0:250)
        write(*,*)'  please enter v, x and y ( z=x+iy )'
        read(*,*)v,x,y
        z=cmplx(x,y)
        n=int(v)
        v0=v-n
        write(*,25)v0,x,y
        if (n.le.8) then
           ns=1
        else
           write(*,*)'  please enter order step ns'
           read(*,*)ns
        endif
        call cikvb(v,z,vm,cbi,cdi,cbk,cdk)
        nm=int(vm)
        write(*,*)
        write(*,*)'   n      re[iv(z)]       im[iv(z)] ',
     &            '     re[iv''(z)]      im[iv''(z)] '
        write(*,*)' ---------------------------------------',
     &            '------------------------------'
        do 10 k=0,nm,ns
10         write(*,20) k,cbi(k),cdi(k)
        write(*,*)
        write(*,*)'   n      re[kv(z)]       im[kv(z)] ',
     &            '     re[kv''(z)]      im[kv''(z)] '
        write(*,*)' ---------------------------------------',
     &            '------------------------------'
        do 15 k=0,nm,ns
15         write(*,20)k,cbk(k),cdk(k)
20      format(1x,i4,1x,4d16.8)
25      format(8x,'v= n+v0',',   ','v0 =',f6.2,',   ','z =',f5.1,
     &        '+ i',f5.1)
        end


        subroutine cikvb(v,z,vm,cbi,cdi,cbk,cdk)
c
c       ===========================================================
c       purpose: compute the modified bessel functions iv(z), kv(z)
c                and their derivatives for an arbitrary order and
c                complex argument
c       input :  z --- complex argument z
c                v --- real order of iv(z) and kv(z)
c                      ( v =n+v0, n = 0,1,2,..., 0 ף v0 < 1 )
c       output:  cbi(n) --- in+v0(z)
c                cdi(n) --- in+v0'(z)
c                cbk(n) --- kn+v0(z)
c                cdk(n) --- kn+v0'(z)
c                vm --- highest order computed
c       routines called:
c            (1) gamma for computing the gamma function
c            (2) msta1 and msta2 for computing the starting
c                point for backward recurrence
c       ===========================================================
c
        implicit double precision (a,d-h,o-y)
        implicit complex*16 (c,z)
        dimension cbi(0:*),cdi(0:*),cbk(0:*),cdk(0:*)
        z1=z
        z2=z*z
        a0=cdabs(z)
        pi=3.141592653589793d0
        ci=(0.0d0,1.0d0)
        n=int(v)
        v0=v-n
        piv=pi*v0
        vt=4.0d0*v0*v0
        if (n.eq.0) n=1
        if (a0.lt.1.0d-100) then
           do 10 k=0,n
              cbi(k)=0.0d0
              cdi(k)=0.0d0
              cbk(k)=-1.0d+300
10            cdk(k)=1.0d+300
           if (v0.eq.0.0) then
              cbi(0)=(1.0d0,0.0d0)
              cdi(1)=(0.5d0,0.0d0)
           endif
           vm=v
           return
        endif
        k0=14
        if (a0.ge.35.0) k0=10
        if (a0.ge.50.0) k0=8
        if (real(z).lt.0.0) z1=-z
        if (a0.lt.18.0) then
           if (v0.eq.0.0) then
              ca1=(1.0d0,0.0d0)
           else
              v0p=1.0d0+v0
              call gamma(v0p,gap)
              ca1=(0.5d0*z1)**v0/gap
           endif
           ci0=(1.0d0,0.0d0)
           cr=(1.0d0,0.0d0)
           do 15 k=1,50
              cr=0.25d0*cr*z2/(k*(k+v0))
              ci0=ci0+cr
              if (cdabs(cr/ci0).lt.1.0d-15) go to 20
15         continue
20         cbi0=ci0*ca1
        else
           ca=cdexp(z1)/cdsqrt(2.0d0*pi*z1)
           cs=(1.0d0,0.0d0)
           cr=(1.0d0,0.0d0)
           do 25 k=1,k0
              cr=-0.125d0*cr*(vt-(2.0d0*k-1.0d0)**2.0)/(k*z1)
25            cs=cs+cr
           cbi0=ca*cs
        endif
        m=msta1(a0,200)
        if (m.lt.n) then
           n=m
        else
           m=msta2(a0,n,15)
        endif
        cf2=(0.0d0,0.0d0)
        cf1=(1.0d-100,0.0d0)
        do 30 k=m,0,-1
           cf=2.0d0*(v0+k+1.0d0)/z1*cf1+cf2
           if (k.le.n) cbi(k)=cf
           cf2=cf1
30         cf1=cf
        cs=cbi0/cf
        do 35 k=0,n
35         cbi(k)=cs*cbi(k)
        if (a0.le.9.0) then
           if (v0.eq.0.0) then
              ct=-cdlog(0.5d0*z1)-0.5772156649015329d0
              cs=(0.0d0,0.0d0)
              w0=0.0d0
              cr=(1.0d0,0.0d0)
              do 40 k=1,50
                 w0=w0+1.0d0/k
                 cr=0.25d0*cr/(k*k)*z2
                 cp=cr*(w0+ct)
                 cs=cs+cp
                 if (k.ge.10.and.cdabs(cp/cs).lt.1.0d-15) go to 45
40            continue
45            cbk0=ct+cs
           else
              v0n=1.0d0-v0
              call gamma(v0n,gan)
              ca2=1.0d0/(gan*(0.5d0*z1)**v0)
              ca1=(0.5d0*z1)**v0/gap
              csu=ca2-ca1
              cr1=(1.0d0,0.0d0)
              cr2=(1.0d0,0.0d0)
              do 50 k=1,50
                 cr1=0.25d0*cr1*z2/(k*(k-v0))
                 cr2=0.25d0*cr2*z2/(k*(k+v0))
                 cp=ca2*cr1-ca1*cr2
                 csu=csu+cp
                 if (k.ge.10.and.cdabs(cp/csu).lt.1.0d-15) go to 55
50            continue
55            cbk0=0.5d0*pi*csu/dsin(piv)
           endif
        else
           cb=cdexp(-z1)*cdsqrt(0.5d0*pi/z1)
           cs=(1.0d0,0.0d0)
           cr=(1.0d0,0.0d0)
           do 60 k=1,k0
              cr=0.125d0*cr*(vt-(2.0d0*k-1.0d0)**2.0)/(k*z1)
60            cs=cs+cr
           cbk0=cb*cs
        endif
        cbk(0)=cbk0
        if (real(z).lt.0.0) then
           do 65 k=0,n
              cvk=cdexp((k+v0)*pi*ci)
              if (dimag(z).lt.0.0d0) then
                 cbk(k)=cvk*cbk(k)+pi*ci*cbi(k)
                 cbi(k)=cbi(k)/cvk
              else if (dimag(z).gt.0.0) then
                 cbk(k)=cbk(k)/cvk-pi*ci*cbi(k)
                 cbi(k)=cvk*cbi(k)
              endif
65         continue
        endif
        do 70 k=1,n
           ckk=(1.0d0/z-cbi(k)*cbk(k-1))/cbi(k-1)
           cbk(k)=ckk
70      continue
        cdi(0)=v0/z*cbi(0)+cbi(1)
        cdk(0)=v0/z*cbk(0)-cbk(1)
        do 80 k=1,n
           cdi(k)=-(k+v0)/z*cbi(k)+cbi(k-1)
80         cdk(k)=-(k+v0)/z*cbk(k)-cbk(k-1)
        vm=n+v0
        return
        end


        subroutine gamma(x,ga)
c
c       ==================================================
c       purpose: compute gamma function ג(x)
c       input :  x  --- argument of ג(x)
c                       ( x is not equal to 0,-1,-2,תתת)
c       output:  ga --- ג(x)
c       ==================================================
c
        implicit double precision (a-h,o-z)
        dimension g(26)
        pi=3.141592653589793d0
        if (x.eq.int(x)) then
           if (x.gt.0.0d0) then
              ga=1.0d0
              m1=x-1
              do 10 k=2,m1
10               ga=ga*k
           else
              ga=1.0d+300
           endif
        else
           if (dabs(x).gt.1.0d0) then
              z=dabs(x)
              m=int(z)
              r=1.0d0
              do 15 k=1,m
15               r=r*(z-k)
              z=z-m
           else
              z=x
           endif
           data g/1.0d0,0.5772156649015329d0,
     &          -0.6558780715202538d0, -0.420026350340952d-1,
     &          0.1665386113822915d0,-.421977345555443d-1,
     &          -.96219715278770d-2, .72189432466630d-2,
     &          -.11651675918591d-2, -.2152416741149d-3,
     &          .1280502823882d-3, -.201348547807d-4,
     &          -.12504934821d-5, .11330272320d-5,
     &          -.2056338417d-6, .61160950d-8,
     &          .50020075d-8, -.11812746d-8,
     &          .1043427d-9, .77823d-11,
     &          -.36968d-11, .51d-12,
     &          -.206d-13, -.54d-14, .14d-14, .1d-15/
           gr=g(26)
           do 20 k=25,1,-1
20            gr=gr*z+g(k)
           ga=1.0d0/(gr*z)
           if (dabs(x).gt.1.0d0) then
              ga=ga*r
              if (x.lt.0.0d0) ga=-pi/(x*ga*dsin(pi*x))
           endif
        endif
        return
        end


        integer function msta1(x,mp)
c
c       ===================================================
c       purpose: determine the starting point for backward  
c                recurrence such that the magnitude of    
c                jn(x) at that point is about 10^(-mp)
c       input :  x     --- argument of jn(x)
c                mp    --- value of magnitude
c       output:  msta1 --- starting point   
c       ===================================================
c
        implicit double precision (a-h,o-z)
        a0=dabs(x)
        n0=int(1.1*a0)+1
        f0=envj(n0,a0)-mp
        n1=n0+5
        f1=envj(n1,a0)-mp
        do 10 it=1,20             
           nn=n1-(n1-n0)/(1.0d0-f0/f1)                  
           f=envj(nn,a0)-mp
           if(abs(nn-n1).lt.1) go to 20
           n0=n1
           f0=f1
           n1=nn
 10        f1=f
 20     msta1=nn
        return
        end


        integer function msta2(x,n,mp)
c
c       ===================================================
c       purpose: determine the starting point for backward
c                recurrence such that all jn(x) has mp
c                significant digits
c       input :  x  --- argument of jn(x)
c                n  --- order of jn(x)
c                mp --- significant digit
c       output:  msta2 --- starting point
c       ===================================================
c
        implicit double precision (a-h,o-z)
        a0=dabs(x)
        hmp=0.5d0*mp
        ejn=envj(n,a0)
        if (ejn.le.hmp) then
           obj=mp
           n0=int(1.1*a0)
        else
           obj=hmp+ejn
           n0=n
        endif
        f0=envj(n0,a0)-obj
        n1=n0+5
        f1=envj(n1,a0)-obj
        do 10 it=1,20
           nn=n1-(n1-n0)/(1.0d0-f0/f1)
           f=envj(nn,a0)-obj
           if (abs(nn-n1).lt.1) go to 20
           n0=n1
           f0=f1
           n1=nn
10         f1=f
20      msta2=nn+10
        return
        end

        real*8 function envj(n,x)
        double precision x
        envj=0.5d0*dlog10(6.28d0*n)-n*dlog10(1.36d0*x/n)
        return
        end
