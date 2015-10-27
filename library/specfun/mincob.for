        program mincob
c
c       =========================================================
c       purpose: this program computes the incomplete beta
c                function ix(a,b) using subroutine incob
c       input :  a --- parameter
c                b --- parameter
c                x --- argument ( 0 ף x ף 1 )
c       output:  bix --- ix(a,b)
c       example:
c                  a       b       x       ix(a,b)
c                -----------------------------------
c                 1.0     3.0     .25     .57812500
c       =========================================================
c
        implicit double precision (a-h,o-z)
        write(*,*)'please enter a, b and x ( 0 ף x ף 1 )'
        read(*,*)a,b,x
        write(*,*)
        write(*,*)'   a       b       x       ix(a,b)'
        write(*,*)' -----------------------------------'
        call incob(a,b,x,bix)
        write(*,10)a,b,x,bix
10      format(1x,f5.1,3x,f5.1,3x,f5.2,f14.8)
        end


        subroutine incob(a,b,x,bix)
c
c       ========================================================
c       purpose: compute the incomplete beta function ix(a,b)
c       input :  a --- parameter
c                b --- parameter
c                x --- argument ( 0 ף x ף 1 )
c       output:  bix --- ix(a,b)
c       routine called: beta for computing beta function b(p,q)
c       ========================================================
c
        implicit double precision (a-h,o-z)
        dimension dk(51),fk(51)
        s0=(a+1.0d0)/(a+b+2.0d0)
        call beta(a,b,bt)
        if (x.le.s0) then
           do 10 k=1,20
10            dk(2*k)=k*(b-k)*x/(a+2.0d0*k-1.0d0)/(a+2.0d0*k)
           do 15 k=0,20
15            dk(2*k+1)=-(a+k)*(a+b+k)*x/(a+2.d0*k)/(a+2.0*k+1.0)
           t1=0.0d0
           do 20 k=20,1,-1
20            t1=dk(k)/(1.0d0+t1)
           ta=1.0d0/(1.0d0+t1)
           bix=x**a*(1.0d0-x)**b/(a*bt)*ta
        else
           do 25 k=1,20
25            fk(2*k)=k*(a-k)*(1.0d0-x)/(b+2.*k-1.0)/(b+2.0*k)
           do 30 k=0,20
30            fk(2*k+1)=-(b+k)*(a+b+k)*(1.d0-x)/
     &                   (b+2.d0*k)/(b+2.d0*k+1.d0)
           t2=0.0d0
           do 35 k=20,1,-1
35            t2=fk(k)/(1.0d0+t2)
           tb=1.0d0/(1.0d0+t2)
           bix=1.0d0-x**a*(1.0d0-x)**b/(b*bt)*tb
        endif
        return
        end


        subroutine beta(p,q,bt)
c
c       ==========================================
c       purpose: compute the beta function b(p,q)
c       input :  p --- parameter  ( p > 0 )
c                q --- parameter  ( q > 0 )
c       output:  bt --- b(p,q)
c       routine called: gamma for computing ג(x)
c       ==========================================
c
        implicit double precision (a-h,o-z)
        call gamma(p,gp)
        call gamma(q,gq)
        ppq=p+q
        call gamma(ppq,gpq)
        bt=gp*gq/gpq
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
