        program mlgama
c
c       ===================================================
c       purpose: this program computes the gamma function
c                â(x) for x > 0 using subroutine lgama
c       examples:
c                  x           â(x)
c                -------------------------
c                 0.5     .1772453851d+01
c                 2.5     .1329340388d+01
c                 5.0     .2400000000d+02
c                 7.5     .1871254306d+04
c                10.0     .3628800000d+06
c       ===================================================
c
        implicit double precision (g,x)
        write(*,*)'   x           â(x)'
        write(*,*)' -------------------------'
        do 10 l=0,20,5
           x=0.5d0*l
           if (l.eq.0) x=0.5
           call lgama(1,x,gl)
           write(*,20)x,gl
10      continue
        write(*,*) 'please enter x:'
        read(*,*) x
        call lgama(1,x,gl)
        write(*,20)x,gl
20      format(1x,f5.1,d20.10)
        end


        subroutine lgama(kf,x,gl)
c
c       ==================================================
c       purpose: compute gamma function â(x) or ln[â(x)]
c       input:   x  --- argument of â(x) ( x > 0 )
c                kf --- function code
c                       kf=1 for â(x); kf=0 for ln[â(x)]
c       output:  gl --- â(x) or ln[â(x)]
c       ==================================================
c
        implicit double precision (a-h,o-z)
        dimension a(10)
        data a/8.333333333333333d-02,-2.777777777777778d-03,
     &         7.936507936507937d-04,-5.952380952380952d-04,
     &         8.417508417508418d-04,-1.917526917526918d-03,
     &         6.410256410256410d-03,-2.955065359477124d-02,
     &         1.796443723688307d-01,-1.39243221690590d+00/
        x0=x
        if (x.eq.1.0.or.x.eq.2.0) then
           gl=0.0d0
           go to 20
        else if (x.le.7.0) then
           n=int(7-x)
           x0=x+n
        endif
        x2=1.0d0/(x0*x0)
        xp=6.283185307179586477d0
        gl0=a(10)
        do 10 k=9,1,-1
10         gl0=gl0*x2+a(k)
        gl=gl0/x0+0.5d0*dlog(xp)+(x0-.5d0)*dlog(x0)-x0
        if (x.le.7.0) then
           do 15 k=1,n
              gl=gl-dlog(x0-1.0d0)
15            x0=x0-1.0d0
        endif
20      if (kf.eq.1) gl=dexp(gl)
        return
        end
