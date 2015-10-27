        program mpbdv
c
c       =========================================================
c       purpose: this program computes the parabolic cylinder 
c                functions dv(x) and their derivatives using
c                subroutine pbdv
c       input:   x --- argument of dv(x)
c                v --- order of dv(x)
c       output:  dv(na) --- dn+v0(x)
c                dp(na) --- dn+v0'(x)
c                ( na = |n|, n = int(v), v0 = v-n, |v0| < 1
c                  n = 0,ס1,ס2,תתת, |n| ף 100 )
c                pdf --- dv(x)
c                pdd --- dv'(x)
c       example: v = 5.5,  x =10.0,  v0 = 0.5,  n = 0,1,...,5
c
c                  n+v0      dv(x)           dv'(x)
c                ---------------------------------------
c                  0.5   .43971930d-10  -.21767183d-09
c                  1.5   .43753148d-09  -.21216995d-08
c                  2.5   .43093569d-08  -.20452956d-07
c                  3.5   .41999741d-07  -.19491595d-06
c                  4.5   .40491466d-06  -.18355745d-05
c                  5.5   .38601477d-05  -.17073708d-04
c
c                dv(x)= .38601477d-05,  dv'(x)=-.17073708d-04
c       =========================================================
c
        implicit double precision (a-h,o-z)
        dimension dv(0:100),dp(0:100)
        write(*,*)'please enter v and  x '
        read(*,*)v,x
        write(*,20)v,x
        nv=int(v)
        v0=v-nv
        na=abs(nv)
        call pbdv(v,x,dv,dp,pdf,pdd)
        write(*,*)
        write(*,*)'   v       dv(x)           dv''(x)'
        write(*,*)'---------------------------------------'
        do 10 k=0,na
           vk=k*isign(1,nv)+v0
10         write(*,30)vk,dv(k),dp(k)
        write(*,*)
        write(*,40)v,pdf,pdd
20      format(1x,'v =',f6.2,',    ','x =',f6.2)
30      format(1x,f5.1,2d16.8)
40      format(1x,'v =',f5.1,',  dv(x)=',d14.8,',   dv''(x)=',d14.8)
        end


        subroutine pbdv(v,x,dv,dp,pdf,pdd)
c
c       ====================================================
c       purpose: compute parabolic cylinder functions dv(x)
c                and their derivatives
c       input:   x --- argument of dv(x)
c                v --- order of dv(x)
c       output:  dv(na) --- dn+v0(x)
c                dp(na) --- dn+v0'(x)
c                ( na = |n|, v0 = v-n, |v0| < 1, 
c                  n = 0,ס1,ס2,תתת )
c                pdf --- dv(x)
c                pdd --- dv'(x)
c       routines called:
c             (1) dvsa for computing dv(x) for small |x|
c             (2) dvla for computing dv(x) for large |x|
c       ====================================================
c
        implicit double precision (a-h,o-z)
        dimension dv(0:*),dp(0:*)
        xa=dabs(x)
        vh=v
        v=v+dsign(1.0d0,v)
        nv=int(v)
        v0=v-nv
        na=abs(nv)
        ep=dexp(-.25d0*x*x)
        if (na.ge.1) ja=1
        if (v.ge.0.0) then
           if (v0.eq.0.0) then
              pd0=ep
              pd1=x*ep
           else
              do 10 l=0,ja
                 v1=v0+l
                 if (xa.le.5.8) call dvsa(v1,x,pd1)
                 if (xa.gt.5.8) call dvla(v1,x,pd1)
                 if (l.eq.0) pd0=pd1
10            continue
           endif
           dv(0)=pd0
           dv(1)=pd1
           do 15 k=2,na
              pdf=x*pd1-(k+v0-1.0d0)*pd0
              dv(k)=pdf
              pd0=pd1
15            pd1=pdf
        else
           if (x.le.0.0) then
              if (xa.le.5.8d0)  then
                 call dvsa(v0,x,pd0)
                 v1=v0-1.0d0
                 call dvsa(v1,x,pd1)
              else
                 call dvla(v0,x,pd0)
                 v1=v0-1.0d0
                 call dvla(v1,x,pd1)
              endif
              dv(0)=pd0
              dv(1)=pd1
              do 20 k=2,na
                 pd=(-x*pd1+pd0)/(k-1.0d0-v0)
                 dv(k)=pd
                 pd0=pd1
20               pd1=pd
           else if (x.le.2.0) then
              v2=nv+v0
              if (nv.eq.0) v2=v2-1.0d0
              nk=int(-v2)
              call dvsa(v2,x,f1)
              v1=v2+1.0d0
              call dvsa(v1,x,f0)
              dv(nk)=f1
              dv(nk-1)=f0
              do 25 k=nk-2,0,-1
                 f=x*f0+(k-v0+1.0d0)*f1
                 dv(k)=f
                 f1=f0
25               f0=f
           else
              if (xa.le.5.8) call dvsa(v0,x,pd0)
              if (xa.gt.5.8) call dvla(v0,x,pd0)
              dv(0)=pd0
              m=100+na
              f1=0.0d0
              f0=1.0d-30
              do 30 k=m,0,-1
                 f=x*f0+(k-v0+1.0d0)*f1
                 if (k.le.na) dv(k)=f
                 f1=f0
30               f0=f
              s0=pd0/f
              do 35 k=0,na
35               dv(k)=s0*dv(k)
           endif
        endif
        do 40 k=0,na-1
           v1=abs(v0)+k
           if (v.ge.0.0d0) then
              dp(k)=0.5d0*x*dv(k)-dv(k+1)
           else
              dp(k)=-0.5d0*x*dv(k)-v1*dv(k+1)
           endif
40      continue
        pdf=dv(na-1)
        pdd=dp(na-1)
        v=vh
        return
        end


        subroutine dvsa(va,x,pd)
c
c       ===================================================
c       purpose: compute parabolic cylinder function dv(x)
c                for small argument
c       input:   x  --- argument
c                va --- order
c       output:  pd --- dv(x)
c       routine called: gamma for computing ג(x)
c       ===================================================
c
        implicit double precision (a-h,o-z)
        eps=1.0d-15
        pi=3.141592653589793d0
        sq2=dsqrt(2.0d0)
        ep=dexp(-.25d0*x*x)
        va0=0.5d0*(1.0d0-va)
        if (va.eq.0.0) then
           pd=ep
        else
           if (x.eq.0.0) then
              if (va0.le.0.0.and.va0.eq.int(va0)) then
                 pd=0.0d0
              else
                 call gamma(va0,ga0)
                 pd=dsqrt(pi)/(2.0d0**(-.5d0*va)*ga0)
              endif
           else
              call gamma(-va,g1)
              a0=2.0d0**(-0.5d0*va-1.0d0)*ep/g1
              vt=-.5d0*va
              call gamma(vt,g0)
              pd=g0
              r=1.0d0
              do 10 m=1,250
                 vm=.5d0*(m-va)
                 call gamma(vm,gm)
                 r=-r*sq2*x/m
                 r1=gm*r
                 pd=pd+r1
                 if (dabs(r1).lt.dabs(pd)*eps) go to 15
10            continue
15            pd=a0*pd
           endif
        endif
        return
        end


        subroutine dvla(va,x,pd)
c
c       ====================================================
c       purpose: compute parabolic cylinder functions dv(x)
c                for large argument
c       input:   x  --- argument
c                va --- order
c       output:  pd --- dv(x)
c       routines called:
c             (1) vvla for computing vv(x) for large |x|
c             (2) gamma for computing ג(x)
c       ====================================================
c
        implicit double precision (a-h,o-z)
        pi=3.141592653589793d0           
        eps=1.0d-12
        ep=dexp(-.25*x*x)
        a0=dabs(x)**va*ep
        r=1.0d0
        pd=1.0d0
        do 10 k=1,16
           r=-0.5d0*r*(2.0*k-va-1.0)*(2.0*k-va-2.0)/(k*x*x)
           pd=pd+r
           if (dabs(r/pd).lt.eps) go to 15
10      continue
15      pd=a0*pd
        if (x.lt.0.0d0) then
            x1=-x
            call vvla(va,x1,vl)
            call gamma(-va,gl)
            pd=pi*vl/gl+dcos(pi*va)*pd
        endif
        return
        end


        subroutine vvla(va,x,pv)
c
c       ===================================================
c       purpose: compute parabolic cylinder function vv(x)
c                for large argument
c       input:   x  --- argument
c                va --- order
c       output:  pv --- vv(x)
c       routines called:
c             (1) dvla for computing dv(x) for large |x|
c             (2) gamma for computing ג(x)
c       ===================================================
c
        implicit double precision (a-h,o-z)
        pi=3.141592653589793d0
        eps=1.0d-12
        qe=dexp(0.25*x*x)
        a0=dabs(x)**(-va-1.0d0)*dsqrt(2.0d0/pi)*qe
        r=1.0d0
        pv=1.0d0
        do 10 k=1,18
           r=0.5d0*r*(2.0*k+va-1.0)*(2.0*k+va)/(k*x*x)
           pv=pv+r
           if (dabs(r/pv).lt.eps) go to 15
10      continue
15      pv=a0*pv
        if (x.lt.0.0d0) then
           x1=-x
           call dvla(va,x1,pdl)
           call gamma(-va,gl)
           dsl=dsin(pi*va)*dsin(pi*va)
           pv=dsl*gl/pi*pdl-dcos(pi*va)*pv
        endif
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
