        program mlpmv
c
c       =========================================================
c       purpose: this program computes the associated legendre 
c                function pmv(x) with an integer order and an
c                arbitrary nonnegative degree using subroutine 
c                lpmv
c       input :  x   --- argument of pm(x)  ( -1 ó x ó 1 )
c                m   --- order of pmv(x)
c                v   --- degree of pmv(x)
c       output:  pmv --- pmv(x)
c       example:    m = 4,  x = 0.5
c                    v          pmv(x)
c                 -----------------------
c                   1.5       .46218726
c                   1.6       .48103143
c                   1.7       .45031429
c                   1.8       .36216902
c                   1.9       .21206446
c                   2.0       .00000000
c                   2.5     -1.51996235
c       =========================================================
c
        implicit double precision (p,v,x)
        write(*,*)'please enter m,v,x = ?'
        read(*,*) m,v,x
        write(*,20)m,x
        write(*,*)
        write(*,*)'     v        pmv(x)'
        write(*,*)'  -----------------------'
        call lpmv(v,m,x,pmv)
        write(*,10)v,pmv
10      format(3x,f5.1,e16.8)
20      format(3x,'m =',i2,',    ','x =',f6.2)
        end


        subroutine lpmv(v,m,x,pmv)
c
c       =======================================================
c       purpose: compute the associated legendre function
c                pmv(x) with an integer order and an arbitrary 
c                nonnegative degree v
c       input :  x   --- argument of pm(x)  ( -1 ó x ó 1 )
c                m   --- order of pmv(x)
c                v   --- degree of pmv(x)
c       output:  pmv --- pmv(x)
c       routine called:  psi for computing psi function
c       =======================================================
c
        implicit double precision (a-h,o-z)
        pi=3.141592653589793d0
        el=.5772156649015329d0
        eps=1.0d-14
        nv=int(v)
        v0=v-nv
        if (x.eq.-1.0d0.and.v.ne.nv) then
           if (m.eq.0) pmv=-1.0d+300
           if (m.ne.0) pmv=1.0d+300
           return
        endif
        c0=1.0d0
        if (m.ne.0) then
           rg=v*(v+m)
           do 10 j=1,m-1
10            rg=rg*(v*v-j*j)
           xq=dsqrt(1.0d0-x*x)
           r0=1.0d0
           do 15 j=1,m
15            r0=.5d0*r0*xq/j
           c0=r0*rg
        endif
        if (v0.eq.0.0d0) then
           pmv=1.0d0
           r=1.0d0
           do 20 k=1,nv-m
              r=0.5d0*r*(-nv+m+k-1.0d0)*(nv+m+k)/(k*(k+m))
     &          *(1.0d0+x)
20            pmv=pmv+r
           pmv=(-1)**nv*c0*pmv
        else
           if (x.ge.-0.35d0) then
              pmv=1.0d0
              r=1.0d0
              do 25 k=1,100
                 r=0.5d0*r*(-v+m+k-1.0d0)*(v+m+k)/(k*(m+k))*(1.0d0-x)
                 pmv=pmv+r
                 if (k.gt.12.and.dabs(r/pmv).lt.eps) go to 30
25            continue
30            pmv=(-1)**m*c0*pmv
           else
              vs=dsin(v*pi)/pi
              pv0=0.0d0
              if (m.ne.0) then
                 qr=dsqrt((1.0d0-x)/(1.0d0+x))
                 r2=1.0d0
                 do 35 j=1,m
35                  r2=r2*qr*j
                 s0=1.0d0
                 r1=1.0d0
                 do 40 k=1,m-1
                    r1=0.5d0*r1*(-v+k-1)*(v+k)/(k*(k-m))*(1.0d0+x)
40                  s0=s0+r1
                 pv0=-vs*r2/m*s0
              endif
              call psi(v,psv)
              pa=2.0d0*(psv+el)+pi/dtan(pi*v)+1.0d0/v
              s1=0.0d0
              do 45 j=1,m
45               s1=s1+(j*j+v*v)/(j*(j*j-v*v))
              pmv=pa+s1-1.0d0/(m-v)+dlog(0.5d0*(1.0d0+x))
              r=1.0d0
              do 60 k=1,100
                 r=0.5d0*r*(-v+m+k-1.0d0)*(v+m+k)/(k*(k+m))*(1.0d0+x)
                 s=0.0d0
                 do 50 j=1,m
50                  s=s+((k+j)**2+v*v)/((k+j)*((k+j)**2-v*v))
                 s2=0.0d0
                 do 55 j=1,k
55                  s2=s2+1.0d0/(j*(j*j-v*v))
                 pss=pa+s+2.0d0*v*v*s2-1.0d0/(m+k-v)
     &               +dlog(0.5d0*(1.0d0+x))
                 r2=pss*r
                 pmv=pmv+r2
                 if (dabs(r2/pmv).lt.eps) go to 65
60            continue
65            pmv=pv0+pmv*vs*c0
           endif
        endif
        return
        end


        subroutine psi(x,ps)
c
c       ======================================
c       purpose: compute psi function
c       input :  x  --- argument of psi(x)
c       output:  ps --- psi(x)
c       ======================================
c
        implicit double precision (a-h,o-z)
        xa=dabs(x)
        pi=3.141592653589793d0
        el=.5772156649015329d0
        s=0.0d0
        if (x.eq.int(x).and.x.le.0.0) then
           ps=1.0d+300
           return
        else if (xa.eq.int(xa)) then
           n=xa
           do 10 k=1 ,n-1
10            s=s+1.0d0/k
           ps=-el+s
        else if (xa+0.5.eq.int(xa+0.5)) then
           n=xa-.5
           do 20 k=1,n
20            s=s+1.0/(2.0d0*k-1.0d0)
           ps=-el+2.0d0*s-1.386294361119891d0
        else
           if (xa.lt.10.0) then
              n=10-int(xa)
              do 30 k=0,n-1
30               s=s+1.0d0/(xa+k)
              xa=xa+n
           endif
           x2=1.0d0/(xa*xa)
           a1=-.8333333333333d-01
           a2=.83333333333333333d-02
           a3=-.39682539682539683d-02
           a4=.41666666666666667d-02
           a5=-.75757575757575758d-02
           a6=.21092796092796093d-01
           a7=-.83333333333333333d-01
           a8=.4432598039215686d0
           ps=dlog(xa)-.5d0/xa+x2*(((((((a8*x2+a7)*x2+
     &        a6)*x2+a5)*x2+a4)*x2+a3)*x2+a2)*x2+a1)
           ps=ps-s
        endif
        if (x.lt.0.0) ps=ps-pi*dcos(pi*x)/dsin(pi*x)-1.0d0/x
        return
        end
