        program mcisib
c
c       ========================================================
c       purpose: this program computes the cosine and sine 
c                integrals using subroutine cisib
c       input :  x  --- argument of ci(x) and si(x)
c       output:  ci --- ci(x)
c                si --- si(x)
c       example:
c
c                   x        ci(x)           si(x)
c                ------------------------------------
c                  0.0    - ì                 0
c                  5.0    -.190030d+00      1.549931
c                 10.0    -.454563d-01      1.658348
c                 20.0     .444201d-01      1.548241
c                 30.0    -.330326d-01      1.566757
c                 40.0     .190201d-01      1.586985
c       ========================================================
c
        implicit double precision (a-h,o-z)
        write(*,*)'please enter x '
        read(*,*) x
        write(*,*)'   x        ci(x)           si(x)'
        write(*,*)'------------------------------------'
        call cisib(x,ci,si)
        if (x.ne.0.0d0) write(*,10)x,ci,si
        if (x.eq.0.0d0) write(*,20)
10      format(1x,f5.1,d16.6,f14.6)
20      format(3x,' .0',3x,' - ì',17x,'0')
        end


        subroutine cisib(x,ci,si)
c
c       =============================================
c       purpose: compute cosine and sine integrals
c                si(x) and ci(x) ( x ò 0 )
c       input :  x  --- argument of ci(x) and si(x)
c       output:  ci --- ci(x)
c                si --- si(x)
c       =============================================
c
        implicit double precision (a-h,o-z)
        x2=x*x
        if (x.eq.0.0) then
           ci=-1.0d+300
           si=0.0d0
        else if (x.le.1.0d0) then
           ci=((((-3.0d-8*x2+3.10d-6)*x2-2.3148d-4)
     &        *x2+1.041667d-2)*x2-0.25)*x2+0.577215665d0+log(x)
           si=((((3.1d-7*x2-2.834d-5)*x2+1.66667d-003)
     &        *x2-5.555556d-002)*x2+1.0)*x
        else
           fx=((((x2+38.027264d0)*x2+265.187033d0)*x2
     &        +335.67732d0)*x2+38.102495d0)/((((x2
     &        +40.021433d0)*x2+322.624911d0)*x2
     &        +570.23628d0)*x2+157.105423d0)
           gx=((((x2+42.242855d0)*x2+302.757865d0)*x2
     &        +352.018498d0)*x2+21.821899d0)/((((x2
     &        +48.196927d0)*x2+482.485984d0)*x2
     &        +1114.978885d0)*x2+449.690326d0)/x
           ci=fx*sin(x)/x-gx*cos(x)/x
           si=1.570796327d0-fx*cos(x)/x-gx*sin(x)/x
        endif
        return
        end
