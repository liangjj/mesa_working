        program msphy
c
c       ========================================================
c       purpose: this program computes the spherical bessel 
c                functions yn(x) and yn'(x) using subroutine
c                sphy
c       input :  x --- argument of yn(x) ( x ע 0 )
c                n --- order of yn(x) ( n = 0,1,תתת, ף 250 )
c       output:  sy(n) --- yn(x)
c                dy(n) --- yn'(x)
c       example:   x = 10.0
c                  n          yn(x)               yn'(x)
c                --------------------------------------------
c                  0     .8390715291d-01    -.6279282638d-01
c                  1     .6279282638d-01     .7134858763d-01
c                  2    -.6506930499d-01     .8231361788d-01
c                  3    -.9532747888d-01    -.2693831344d-01
c                  4    -.1659930220d-02    -.9449751377d-01
c                  5     .9383354168d-01    -.5796005523d-01
c       ========================================================
c
        implicit double precision (a-h,o-z)
        dimension sy(0:250),dy(0:250)
        write(*,*)'please enter n and x '
        read(*,*)n,x
        write(*,30)n,x
        if (n.le.10) then
           ns=1
        else
           write(*,*)'please enter order step ns'
           read(*,*)ns
        endif
        call sphy(n,x,nm,sy,dy)
        write(*,*)
        write(*,*)'  n          yn(x)               yn''(x)'
        write(*,*)'--------------------------------------------'
        do 10 k=0,nm,ns
10         write(*,20)k,sy(k),dy(k)
20      format(1x,i3,2d20.10)
30      format(3x,6hnmax =,i3,',     ',2hx=,f6.1)
        end


        subroutine sphy(n,x,nm,sy,dy)
c
c       ======================================================
c       purpose: compute spherical bessel functions yn(x) and
c                their derivatives
c       input :  x --- argument of yn(x) ( x ע 0 )
c                n --- order of yn(x) ( n = 0,1,תתת )
c       output:  sy(n) --- yn(x)
c                dy(n) --- yn'(x)
c                nm --- highest order computed
c       ======================================================
c
        implicit double precision (a-h,o-z)
        dimension sy(0:n),dy(0:n)
        nm=n
        if (x.lt.1.0d-60) then
           do 10 k=0,n
              sy(k)=-1.0d+300
10            dy(k)=1.0d+300
           return
        endif
        sy(0)=-dcos(x)/x
        sy(1)=(sy(0)-dsin(x))/x
        f0=sy(0)
        f1=sy(1)
        do 15 k=2,n
           f=(2.0d0*k-1.0d0)*f1/x-f0
           sy(k)=f
           if (dabs(f).ge.1.0d+300) go to 20              
           f0=f1
15         f1=f
20      nm=k-1
           dy(0)=(dsin(x)+dcos(x)/x)/x
           do 25 k=1,nm
25            dy(k)=sy(k-1)-(k+1.0d0)*sy(k)/x
        return
        end
