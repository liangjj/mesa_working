        program mstvlv
c
c       ======================================================
c       purpose:  this program computes the modified struve 
c                 function lv(x) for an arbitrary order v
c                 using subroutine stvlv
c       input :   v   --- order of lv(x)  ( |v| ף 20 )
c                 x   --- argument of lv(x) ( x ע 0 )
c       output:   slv --- lv(x)
c       example:  x = 10.0
c                   v          lv(x)
c                 ------------------------
c                  0.5     .27785323d+04
c                  1.5     .24996698d+04
c                  2.5     .20254774d+04
c                  3.5     .14816746d+04
c                  4.5     .98173460d+03
c                  5.5     .59154277d+03
c       =====================================================
c
        implicit double precision (a-h,o-z)
        write(*,*)'please enter v and x '
        read(*,*)v,x
        write(*,30)v,x
        write(*,*)
        write(*,*)'   v          lv(x)'
        write(*,*)' ------------------------'
        call stvlv(v,x,slv)
        write(*,20)v,slv
20      format(1x,f5.1,d18.8)
30      format(1x,'v =',f5.1,6x,'x =',f5.1)
        end


        subroutine stvlv(v,x,slv)
c
c       ======================================================
c       purpose:  compute modified struve function lv(x) with
c                 an arbitrary order v
c       input :   v   --- order of lv(x)  ( |v| ף 20 )
c                 x   --- argument of lv(x) ( x ע 0 )
c       output:   slv --- lv(x)
c       routine called: gamma to compute the gamma function
c       ======================================================
c
        implicit double precision (a-h,o-z)
        pi=3.141592653589793d0
        if (x.eq.0.0d0) then
           if (v.gt.-1.0.or.int(v)-v.eq.0.5d0) then
              slv=0.0d0
           else if (v.lt.-1.0d0) then
              slv=(-1)**(int(0.5d0-v)-1)*1.0d+300
           else if (v.eq.-1.0d0) then
              slv=2.0d0/pi
           endif
           return
        endif
        if (x.le.40.0d0) then
           v0=v+1.5d0
           call gamma(v0,ga)
           s=2.0d0/(dsqrt(pi)*ga)
           r1=1.0d0
           do 10 k=1,100
              va=k+1.5d0
              call gamma(va,ga)
              vb=v+k+1.5d0
              call gamma(vb,gb)
              r1=r1*(0.5d0*x)**2
              r2=r1/(ga*gb)
              s=s+r2
              if (dabs(r2/s).lt.1.0d-12) go to 15
10         continue
15         slv=(0.5d0*x)**(v+1.0d0)*s
        else
           sa=-1.0d0/pi*(0.5d0*x)**(v-1.0)
           v0=v+0.5d0
           call gamma(v0,ga)
           s=-dsqrt(pi)/ga
           r1=-1.0d0
           do 20 k=1,12
              va=k+0.5d0
              call gamma(va,ga)
              vb=-k+v+0.5d0
              call gamma(vb,gb)
              r1=-r1/(0.5d0*x)**2
              s=s+r1*ga/gb
20         continue
           s0=sa*s
           u=dabs(v)
           n=int(u)
           u0=u-n
           do 35 l=0,1
              vt=u0+l
              r=1.0d0
              biv=1.0d0
              do 25 k=1,16
                 r=-0.125*r*(4.0*vt*vt-(2.0*k-1.0d0)**2)/(k*x)
                 biv=biv+r
                 if (dabs(r/biv).lt.1.0d-12) go to 30
25            continue
30            if (l.eq.0) biv0=biv
35         continue
           bf0=biv0
           bf1=biv
           do 40 k=2,n
              bf=-2.0d0*(k-1.0+u0)/x*bf1+bf0
              bf0=bf1
40            bf1=bf
           if (n.eq.0) biv=biv0
           if (n.gt.1) biv=bf
           slv=dexp(x)/dsqrt(2.0d0*pi*x)*biv+s0
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
