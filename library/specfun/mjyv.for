        program mjyv
c
c       ============================================================
c       purpose: this program computes bessel functions jv(x) and
c                yv(x) and their derivatives using subroutine jyv
c       input :  x --- argument of jv(x) and yv(x)
c                v --- order of jv(x) and yv(x)
c                      ( v = n+v0,  0 ף n ף 250, 0 ף v0 < 1 )
c       output:  bj(n) --- jn+v0(x)
c                dj(n) --- jn+v0'(x)
c                by(n) --- yn+v0(x)
c                dy(n) --- yn+v0'(x)
c       example: compute jv(x) and yv(x) and their derivatives
c                for v = 0.25(1.0)5.25 and x = 10.0
c                computation results:
c
c                v =  5.25,      x = 10.00
c
c        v        jv(x)         jv'(x)        yv(x)         yv'(x)
c       ------------------------------------------------------------
c        .25   -.20639379    -.13476340     .14493044    -.21381777
c       1.25    .12960355    -.22259423     .21744103     .11775031
c       2.25    .23879467     .07587475    -.09057018     .23781932
c       3.25   -.02214595     .24599211    -.25819761    -.00665596
c       4.25   -.25318954     .08545961    -.07725827    -.22536285
c       5.25   -.19306516    -.15183033     .19252809    -.17833551
c       ============================================================
c
        implicit double precision (a-h,o-z)
        common bj(0:250),dj(0:250),by(0:250),dy(0:250)
        write(*,*)'  please enter v, x '
        read(*,*)v,x
        write(*,20)v,x
        if (v.le.8) then
           ns=1
        else
           write(*,*)'  please enter order step ns'
           read(*,*)ns
        endif
        call jyv(v,x,vm,bj,dj,by,dy)  
        nm=int(vm)
        v0=vm-nm
        write(*,*)
        write(*,*)'   v         jv(x)           jv''(x)',
     &            '          yv(x)           yv''(x)'
        write(*,*)' ---------------------------------------------',
     &            '------------------------'
        do 10 k=0,nm,ns
           vk=k+v0
10         write(*,15)vk,bj(k),dj(k),by(k),dy(k)
15      format(1x,f6.2,4d16.8)
20      format(8x,3hv =,f6.2,',    ',3hx =,f6.2)
        end


        subroutine jyv(v,x,vm,bj,dj,by,dy)
c
c       =======================================================
c       purpose: compute bessel functions jv(x) and yv(x)
c                and their derivatives
c       input :  x --- argument of jv(x) and yv(x)
c                v --- order of jv(x) and yv(x)
c                      ( v = n+v0, 0 ף v0 < 1, n = 0,1,2,... )
c       output:  bj(n) --- jn+v0(x)
c                dj(n) --- jn+v0'(x)
c                by(n) --- yn+v0(x)
c                dy(n) --- yn+v0'(x)
c                vm --- highest order computed
c       routines called:
c            (1) gamma for computing gamma function
c            (2) msta1 and msta2 for computing the starting
c                point for backward recurrence
c       =======================================================
c
        implicit double precision (a-h,o-z)
        dimension bj(0:*),dj(0:*),by(0:*),dy(0:*)
        el=.5772156649015329d0
        pi=3.141592653589793d0
        rp2=.63661977236758d0
        x2=x*x
        n=int(v)
        v0=v-n
        if (x.lt.1.0d-100) then
           do 10 k=0,n
              bj(k)=0.0d0
              dj(k)=0.0d0
              by(k)=-1.0d+300
10            dy(k)=1.0d+300
           if (v0.eq.0.0) then
              bj(0)=1.0d0
              dj(1)=0.5d0
           else
              dj(0)=1.0d+300
           endif
           vm=v  
           return
        endif
        if (x.le.12.0) then
           do 25 l=0,1
              vl=v0+l
              bjvl=1.0d0
              r=1.0d0
              do 15 k=1,40
                 r=-0.25d0*r*x2/(k*(k+vl))
                 bjvl=bjvl+r
                 if (dabs(r).lt.dabs(bjvl)*1.0d-15) go to 20
15            continue
20            vg=1.0d0+vl
              call gamma(vg,ga)
              a=(0.5d0*x)**vl/ga
              if (l.eq.0) bjv0=bjvl*a
              if (l.eq.1) bjv1=bjvl*a
25         continue
        else
           k0=11
           if (x.ge.35.0) k0=10
           if (x.ge.50.0) k0=8
           do 40 j=0,1
              vv=4.0d0*(j+v0)*(j+v0)
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
              xk=x-(0.5d0*(j+v0)+0.25d0)*pi
              a0=dsqrt(rp2/x)
              ck=dcos(xk)
              sk=dsin(xk)
              if (j.eq.0) then
                 bjv0=a0*(px*ck-qx*sk)
                 byv0=a0*(px*sk+qx*ck)
              else if (j.eq.1) then
                 bjv1=a0*(px*ck-qx*sk)
                 byv1=a0*(px*sk+qx*ck)
              endif
40         continue
        endif
        bj(0)=bjv0
        bj(1)=bjv1
        dj(0)=v0/x*bj(0)-bj(1)
        dj(1)=-(1.0d0+v0)/x*bj(1)+bj(0)
        if (n.ge.2.and.n.le.int(0.9*x)) then
           f0=bjv0
           f1=bjv1
           do 45 k=2,n
              f=2.0d0*(k+v0-1.0d0)/x*f1-f0
              bj(k)=f
              f0=f1
45            f1=f
        else if (n.ge.2) then
           m=msta1(x,200)
           if (m.lt.n) then
              n=m
           else
              m=msta2(x,n,15)
           endif
           f2=0.0d0
           f1=1.0d-100
           do 50 k=m,0,-1
              f=2.0d0*(v0+k+1.0d0)/x*f1-f2
              if (k.le.n) bj(k)=f
              f2=f1
50            f1=f
           if (dabs(bjv0).gt.dabs(bjv1)) then
               cs=bjv0/f
           else
               cs=bjv1/f2
           endif
           do 55 k=0,n
55            bj(k)=cs*bj(k)
        endif
        do 60 k=2,n
60         dj(k)=-(k+v0)/x*bj(k)+bj(k-1)
        if (x.le.12.0d0) then
           if (v0.ne.0.0) then
              do 75 l=0,1
                 vl=v0+l
                 bjvl=1.0d0
                 r=1.0d0
                 do 65 k=1,40
                    r=-0.25d0*r*x2/(k*(k-vl))
                    bjvl=bjvl+r
                    if (dabs(r).lt.dabs(bjvl)*1.0d-15) go to 70
65               continue
70               vg=1.0d0-vl
                 call gamma(vg,gb)
                 b=(2.0d0/x)**vl/gb
                 if (l.eq.0) bju0=bjvl*b
                 if (l.eq.1) bju1=bjvl*b
75            continue
              pv0=pi*v0
              pv1=pi*(1.0d0+v0)
              byv0=(bjv0*dcos(pv0)-bju0)/dsin(pv0)
              byv1=(bjv1*dcos(pv1)-bju1)/dsin(pv1)
           else
              ec=dlog(x/2.0d0)+el
              cs0=0.0d0
              w0=0.0d0
              r0=1.0d0
              do 80 k=1,30
                 w0=w0+1.0d0/k
                 r0=-0.25d0*r0/(k*k)*x2
80               cs0=cs0+r0*w0
              byv0=rp2*(ec*bjv0-cs0)
              cs1=1.0d0
              w1=0.0d0
              r1=1.0d0
              do 85 k=1,30
                 w1=w1+1.0d0/k
                 r1=-0.25d0*r1/(k*(k+1))*x2
85               cs1=cs1+r1*(2.0d0*w1+1.0d0/(k+1.0d0))
              byv1=rp2*(ec*bjv1-1.0d0/x-0.25d0*x*cs1)
           endif
        endif
        by(0)=byv0
        by(1)=byv1
        do 90 k=2,n
           byvk=2.0d0*(v0+k-1.0d0)/x*byv1-byv0
           by(k)=byvk
           byv0=byv1
90         byv1=byvk
        dy(0)=v0/x*by(0)-by(1)
        do 95 k=1,n
95         dy(k)=-(k+v0)/x*by(k)+by(k-1)
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
