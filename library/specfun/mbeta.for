        program mbeta
c
c       ====================================================
c       purpose: this program computes the beta function 
c                b(p,q) for p > 0 and q > 0 using 
c                subroutine beta
c       input :  p  --- parameter  ( p > 0 )
c                q  --- parameter  ( q > 0 )
c       output:  bt --- b(p,q)
c       examples:
c                 p       q           b(p,q)
c               ---------------------------------
c                1.5     2.0     .2666666667d+00
c                2.5     2.0     .1142857143d+00
c                1.5     3.0     .1523809524d+00
c       ====================================================
c
        implicit double precision (a-h,o-z)
        write(*,*)'please enter p and q'
        read(*,*)p,q
        write(*,*)
        write(*,*)'    p       q           b(p,q)'
        write(*,*)'  ---------------------------------'
        call beta(p,q,bt)
        write(*,10)p,q,bt
10      format(2x,f5.1,3x,f5.1,d20.10)
        end


        subroutine beta(p,q,bt)
c
c       ==========================================
c       purpose: compute the beta function b(p,q)
c       input :  p  --- parameter  ( p > 0 )
c                q  --- parameter  ( q > 0 )
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
