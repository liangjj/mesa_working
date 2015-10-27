        program mpbvv
c
c       ========================================================
c       purpose: this program computes the parabolic cylinder 
c                functions vv(x) and vv'(x) using subroutine
c                pbvv
c       input:   x --- argument of vv(x)
c                v --- order of vv(x)
c       output:  vv(na) --- vv(x)
c                vp(na) --- vv'(x)
c                ( na = |n|, v = n+v0, n = int(v), |v0| < 1
c                  n = 0,ס1,ס2,תתת, |n| ף 100 )
c                pvf --- vv(x)
c                pvd --- vv'(x)
c       example: v = 5.5,  x =10.0,  v0 = 0.5,  n = 0,1,2,...,5
c
c                  n+v0      vv(x)           vv'(x)
c                ---------------------------------------
c                  0.5   .18522719d+10   .89761157d+10
c                  1.5   .19016268d+09   .90145854d+09
c                  2.5   .19741946d+08   .91452949d+08
c                  3.5   .20733667d+07   .93751130d+07
c                  4.5   .22038231d+06   .97145511d+06
c                  5.5   .23719356d+05   .10178553d+06
c
c                vv(x)= .23719356d+05,  vv'(x)= .10178553d+06
c       =========================================================
c
        implicit double precision (a-h,o-z)
        dimension vv(0:100),vp(0:100)
        write(*,*)'please enter v and  x '
        read(*,*)v,x
        write(*,20)v,x
        nv=int(v)
        v0=v-nv
        na=abs(nv)
        call pbvv(v,x,vv,vp,pvf,pvd)
        write(*,*)
        write(*,*)'   v       vv(x)           vv''(x)'
        write(*,*)'---------------------------------------'
        do 10 k=0,na
           vk=k*isign(1,nv)+v0
10         write(*,30)vk,vv(k),vp(k)
        write(*,*)
        write(*,40)v,pvf,pvd
20      format(1x,'v =',f6.2,',    ','x =',f6.2)
30      format(1x,f5.1,2d16.8)
40      format(1x,'v =',f5.1,',  vv(x)=',d14.8,',   vv''(x)=',d14.8)
        end


        subroutine pbvv(v,x,vv,vp,pvf,pvd)
c
c       ===================================================
c       purpose: compute parabolic cylinder functions vv(x)
c                and their derivatives
c       input:   x --- argument of vv(x)
c                v --- order of vv(x)
c       output:  vv(na) --- vv(x)
c                vp(na) --- vv'(x)
c                ( na = |n|, v = n+v0, |v0| < 1
c                  n = 0,ס1,ס2,תתת )
c                pvf --- vv(x)
c                pvd --- vv'(x)
c       routines called:
c             (1) vvsa for computing vv(x) for small |x|
c             (2) vvla for computing vv(x) for large |x|
c       ===================================================
c
        implicit double precision (a-h,o-z)
        dimension vv(0:*),vp(0:*)
        pi=3.141592653589793d0
        xa=dabs(x)
        vh=v
        v=v+dsign(1.0d0,v)
        nv=int(v)
        v0=v-nv
        na=abs(nv)
        qe=dexp(0.25d0*x*x)
        q2p=dsqrt(2.0d0/pi)
        if (na.ge.1) ja=1
        if (v.le.0.0) then
           if (v0.eq.0.0) then
              if (xa.le.7.5) call vvsa(v0,x,pv0)
              if (xa.gt.7.5) call vvla(v0,x,pv0)
              f0=q2p*qe
              f1=x*f0
              vv(0)=pv0
              vv(1)=f0
              vv(2)=f1
           else
              do 10 l=0,ja
                 v1=v0-l
                 if (xa.le.7.5) call vvsa(v1,x,f1)
                 if (xa.gt.7.5) call vvla(v1,x,f1)
                 if (l.eq.0) f0=f1
10            continue
              vv(0)=f0
              vv(1)=f1
           endif
           kv=2
           if (v0.eq.0.0) kv=3
           do 15 k=kv,na
              f=x*f1+(k-v0-2.0d0)*f0
              vv(k)=f
              f0=f1
15            f1=f
        else
           if (x.ge.0.0.and.x.le.7.5d0) then
              v2=v
              if (v2.lt.1.0) v2=v2+1.0d0
              call vvsa(v2,x,f1)
              v1=v2-1.0d0
              kv=int(v2)
              call vvsa(v1,x,f0)
              vv(kv)=f1
              vv(kv-1)=f0
              do 20 k=kv-2,0,-1
                 f=x*f0-(k+v0+2.0d0)*f1
                 if (k.le.na) vv(k)=f
                 f1=f0
20               f0=f
           else if (x.gt.7.5d0) then
              call vvla(v0,x,pv0)
              m=100+abs(na)
              vv(1)=pv0
              f1=0.0d0
              f0=1.0d-40
              do 25 k=m,0,-1
                 f=x*f0-(k+v0+2.0d0)*f1
                 if (k.le.na) vv(k)=f
                 f1=f0
25               f0=f
              s0=pv0/f
              do 30 k=0,na
30               vv(k)=s0*vv(k)
           else
              if (xa.le.7.5d0) then
                 call vvsa(v0,x,f0)
                 v1=v0+1.0
                 call vvsa(v1,x,f1)
              else
                 call vvla(v0,x,f0)
                 v1=v0+1.0d0
                 call vvla(v1,x,f1)
              endif
              vv(0)=f0
              vv(1)=f1
              do 35 k=2,na
                 f=(x*f1-f0)/(k+v0)
                 vv(k)=f
                 f0=f1
35               f1=f
           endif
        endif
        do 40 k=0,na-1
           v1=v0+k
           if (v.ge.0.0d0) then
              vp(k)=0.5d0*x*vv(k)-(v1+1.0d0)*vv(k+1)
           else
              vp(k)=-0.5d0*x*vv(k)+vv(k+1)
           endif
40      continue
        pvf=vv(na-1)
        pvd=vp(na-1)
        v=vh
        return
        end


        subroutine vvsa(va,x,pv)
c
c       ===================================================
c       purpose: compute parabolic cylinder function vv(x)
c                for small argument
c       input:   x  --- argument
c                va --- order
c       output:  pv --- vv(x)
c       routine called : gamma for computing ג(x)
c       ===================================================
c
        implicit double precision (a-h,o-z)
        eps=1.0d-15
        pi=3.141592653589793d0
        ep=dexp(-.25d0*x*x)
        va0=1.0d0+0.5d0*va
        if (x.eq.0.0) then
           if (va0.le.0.0.and.va0.eq.int(va0).or.va.eq.0.0) then
              pv=0.0d0
           else
              vb0=-0.5d0*va
              sv0=dsin(va0*pi)
              call gamma(va0,ga0)
              pv=2.0d0**vb0*sv0/ga0
           endif
        else
           sq2=dsqrt(2.0d0)
           a0=2.0d0**(-.5d0*va)*ep/(2.0d0*pi)
           sv=dsin(-(va+.5d0)*pi)
           v1=-.5d0*va
           call gamma(v1,g1)
           pv=(sv+1.0d0)*g1
           r=1.0d0
           fac=1.0d0
           do 10 m=1,250
              vm=.5d0*(m-va)
              call gamma(vm,gm)
              r=r*sq2*x/m
              fac=-fac
              gw=fac*sv+1.0d0
              r1=gw*r*gm
              pv=pv+r1
              if (dabs(r1/pv).lt.eps.and.gw.ne.0.0) go to 15
10         continue
15         pv=a0*pv
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
        a0=x**va*ep
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
