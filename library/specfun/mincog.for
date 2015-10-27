        program mincog
c
c       ==========================================================
c       purpose: this program computes the incomplete gamma
c                function r(a,x), ג(a,x) and p(a,x) using
c                subroutine incog
c       input :  a   --- parameter
c                x   --- argument
c       output:  gin --- r(a,x)
c                gim --- ג(a,x)
c                gip --- p(a,x)
c       example:
c            a     x      r(a,x)         ג(a,x)         p(a,x)
c           -------------------------------------------------------
c           3.0   5.0  .17506960d+01  .24930404d+00  .87534798d+00
c       ===========================================================
c

        double precision a,x,gin,gim,gip
        write(*,*)'plese enter a and x'
        read(*,*)a,x
        write(*,*)
        write(*,*)'   a     x      r(a,x)         ג(a,x)         p(a,x)'
        write(*,*)' --------------------------------------------',
     &            '------------'
        call incog(a,x,gin,gim,gip)
        write(*,10)a,x,gin,gim,gip
10      format(1x,f5.1,1x,f5.1,3d15.8)
        end


        subroutine incog(a,x,gin,gim,gip)
c
c       ===================================================
c       purpose: compute the incomplete gamma function
c                r(a,x), ג(a,x) and p(a,x)
c       input :  a   --- parameter ( a ף 170 )
c                x   --- argument 
c       output:  gin --- r(a,x)
c                gim --- ג(a,x)
c                gip --- p(a,x)
c       routine called: gamma for computing ג(x)
c       ===================================================
c
        implicit double precision (a-h,o-z)
        xam=-x+a*dlog(x)
        if (xam.gt.700.0.or.a.gt.170.0) then
           write(*,*)'a and/or x too large'
           stop
        endif
        if (x.eq.0.0) then
           gin=0.0
           call gamma(a,ga)
           gim=ga
           gip=0.0
        else if (x.le.1.0+a) then
           s=1.0d0/a
           r=s
           do 10 k=1,60
              r=r*x/(a+k)
              s=s+r
              if (dabs(r/s).lt.1.0d-15) go to 15
10         continue
15         gin=dexp(xam)*s
           call gamma(a,ga)
           gip=gin/ga
           gim=ga-gin
        else if (x.gt.1.0+a) then
           t0=0.0d0
           do 20 k=60,1,-1
              t0=(k-a)/(1.0d0+k/(x+t0))
20         continue
           gim=dexp(xam)/(x+t0)
           call gamma(a,ga)
           gin=ga-gim
           gip=1.0d0-gim/ga
        endif
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
