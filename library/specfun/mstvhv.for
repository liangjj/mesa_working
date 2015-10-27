        program mstvhv
c
c       ======================================================
c       purpose:  this program computes struve function hv(x) 
c                 for an arbitrary order using subroutine
c                 stvhv
c       input :   v  --- order of hv(x)  ( -8.0 ף v ף 12.5 )
c                 x  --- argument of hv(x) ( x ע 0 )
c       output:   hv --- hv(x)
c       example:  x = 10.0
c                   v           hv(x)
c                 -----------------------
c                   .5     .46402212d+00
c                  1.5     .14452322d+01
c                  2.5     .31234632d+01
c                  3.5     .53730255d+01
c                  4.5     .72083122d+01
c                  5.5     .76851132d+01
c       ======================================================
c
        implicit double precision (a-h,o-z)
        write(*,*)'please v and x '
        read(*,*)v,x
        write(*,30)v,x
        write(*,*)
        write(*,*)'   v           hv(x)'
        write(*,*)' -----------------------'
        call stvhv(v,x,hv)
        write(*,20)v,hv
20      format(1x,f5.1,d18.8)
30      format(1x,'v =',f5.1,6x,'x =',f5.1)
        end


        subroutine stvhv(v,x,hv)
c
c       =====================================================
c       purpose: compute struve function hv(x) with an
c                arbitrary order v
c       input :  v  --- order of hv(x)  ( -8.0 ף v ף 12.5 )
c                x  --- argument of hv(x) ( x ע 0 )
c       output:  hv --- hv(x)
c       routine called: gamma to compute the gamma function
c       =====================================================
c
        implicit double precision (a-h,o-z)
        pi=3.141592653589793d0
        if (x.eq.0.0d0) then
           if (v.gt.-1.0.or.int(v)-v.eq.0.5d0) then
              hv=0.0d0
           else if (v.lt.-1.0d0) then
              hv=(-1)**(int(0.5d0-v)-1)*1.0d+300
           else if (v.eq.-1.0d0) then
              hv=2.0d0/pi
           endif
           return
        endif
        if (x.le.20.0d0) then
           v0=v+1.5d0
           call gamma(v0,ga)
           s=2.0d0/(dsqrt(pi)*ga)
           r1=1.0d0
           do 10 k=1,100
              va=k+1.5d0
              call gamma(va,ga)
              vb=v+k+1.5d0
              call gamma(vb,gb)
              r1=-r1*(0.5d0*x)**2
              r2=r1/(ga*gb)
              s=s+r2
              if (dabs(r2).lt.dabs(s)*1.0d-12) go to 15
10         continue
15         hv=(0.5d0*x)**(v+1.0d0)*s
        else
           sa=(0.5d0*x)**(v-1.0)/pi
           v0=v+0.5d0
           call gamma(v0,ga)
           s=dsqrt(pi)/ga
           r1=1.0d0
           do 20 k=1,12
              va=k+0.5d0
              call gamma(va,ga)
              vb=-k+v+0.5d0
              call gamma(vb,gb)
              r1=r1/(0.5d0*x)**2
              s=s+r1*ga/gb
20         continue
           s0=sa*s
           u=dabs(v)
           n=int(u)
           u0=u-n
           do 35 l=0,1
              vt=4.0d0*(u0+l)**2
              r1=1.0d0
              pu1=1.0d0
              do 25 k=1,12
                 r1=-0.0078125d0*r1*(vt-(4.0*k-3.0d0)**2)*
     &             (vt-(4.0d0*k-1.0)**2)/((2.0d0*k-1.0)*k*x*x)
                 pu1=pu1+r1
25            continue
              qu1=1.0d0
              r2=1.0d0
              do 30 k=1,12
                 r2=-0.0078125d0*r2*(vt-(4.0d0*k-1.0)**2)*
     &             (vt-(4.0d0*k+1.0)**2)/((2.0d0*k+1.0)*k*x*x)
                 qu1=qu1+r2
30            continue
              qu1=0.125d0*(vt-1.0d0)/x*qu1
              if (l.eq.0) then
                 pu0=pu1
                 qu0=qu1
              endif
35         continue
           t0=x-(0.5*u0+0.25d0)*pi
           t1=x-(0.5*u0+0.75d0)*pi
           sr=dsqrt(2.0d0/(pi*x))
           by0=sr*(pu0*dsin(t0)+qu0*dcos(t0))
           by1=sr*(pu1*dsin(t1)+qu1*dcos(t1))
           bf0=by0
           bf1=by1
           do 40 k=2,n
              bf=2.0d0*(k-1.0+u0)/x*bf1-bf0
              bf0=bf1
40            bf1=bf
           if (n.eq.0) byv=by0
           if (n.eq.1) byv=by1
           if (n.gt.1) byv=bf
           hv=byv+s0
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
