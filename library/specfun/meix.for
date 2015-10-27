        program meix
c
c       =========================================================
c       purpose: this program computes the exponential integral 
c                ei(x) using subroutine eix
c       example:
c                  x        ei(x)
c                -----------------------
c                  0    -.10000000+301
c                  1     .18951178e+01
c                  2     .49542344e+01
c                  3     .99338326e+01
c                  4     .19630874e+02
c                  5     .40185275e+02
c       =========================================================
c
        double precision ei,x
        write(*,*)'please enter x '
        read(*,*)x
        write(*,*)
        write(*,*)'   x         ei(x)'
        write(*,*)'------------------------'
        call eix(x,ei)
        write(*,10)x,ei
10      format(1x,f5.1,e18.8)
        end


        subroutine eix(x,ei)
c
c       ============================================
c       purpose: compute exponential integral ei(x)
c       input :  x  --- argument of ei(x)
c       output:  ei --- ei(x) ( x > 0 )
c       ============================================
c
        implicit double precision (a-h,o-z)
        if (x.eq.0.0) then
           ei=-1.0d+300
        else if (x.le.40.0) then
           ei=1.0d0
           r=1.0d0
           do 15 k=1,100
              r=r*k*x/(k+1.0d0)**2
              ei=ei+r
              if (dabs(r/ei).le.1.0d-15) go to 20
15         continue
20         ga=0.5772156649015328d0
           ei=ga+dlog(x)+x*ei
        else
           ei=1.0d0
           r=1.0d0
           do 25 k=1,20
              r=r*k/x
25            ei=ei+r
           ei=dexp(x)/x*ei
        endif
        return
        end
