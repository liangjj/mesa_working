        program msphk
c
c       ======================================================
c       purpose: this program computes the modified spherical 
c                bessel functions kn(x) and kn'(x) using
c                subroutine sphk
c       input :  x --- argument of kn(x)  ( x ò 0 )
c                n --- order of kn(x) ( n ó 250 )
c       output:  sk(n) --- kn(x)
c                dk(n) --- kn'(x)
c       example: x= 10.0
c                  n          kn(x)               kn'(x)
c                --------------------------------------------
c                  0     .7131404291d-05    -.7844544720d-05
c                  1     .7844544720d-05    -.8700313235d-05
c                  2     .9484767707d-05    -.1068997503d-04
c                  3     .1258692857d-04    -.1451953914d-04
c                  4     .1829561771d-04    -.2173473743d-04
c                  5     .2905298451d-04    -.3572740841d-04
c       ======================================================
c
        implicit double precision (a-h,o-z)
        dimension sk(0:250),dk(0:250)
        write(*,*)'please enter n and x '
        read(*,*)n,x
        write(*,30)n,x
        if (n.le.10) then
           ns=1
        else
           write(*,*) 'please enter order step ns'
           read(*,*) ns
        endif
        call sphk(n,x,nm,sk,dk)
        write(*,*)
        write(*,*)'  n          kn(x)               kn''(x)'
        write(*,*)'--------------------------------------------'
        do 10 k=0,nm,ns
10         write(*,20)k,sk(k),dk(k)
20      format(1x,i3,2d20.10)
30      format(3x,'nmax =',i3,',     ','x =',f6.1)
        end


        subroutine sphk(n,x,nm,sk,dk)
c
c       =====================================================
c       purpose: compute modified spherical bessel functions
c                of the second kind, kn(x) and kn'(x)
c       input :  x --- argument of kn(x)  ( x ò 0 )
c                n --- order of kn(x) ( n = 0,1,2,... )
c       output:  sk(n) --- kn(x)
c                dk(n) --- kn'(x)
c                nm --- highest order computed
c       =====================================================
c
        implicit double precision (a-h,o-z)
        dimension sk(0:n),dk(0:n)
        pi=3.141592653589793d0
        nm=n
        if (x.lt.1.0d-60) then
           do 10 k=0,n
              sk(k)=1.0d+300
10            dk(k)=-1.0d+300
           return
        endif
        sk(0)=0.5d0*pi/x*dexp(-x)
        sk(1)=sk(0)*(1.0d0+1.0d0/x)
        f0=sk(0)
        f1=sk(1)
        do 15 k=2,n
           f=(2.0d0*k-1.0d0)*f1/x+f0
           sk(k)=f
           if (dabs(f).gt.1.0d+300) go to 20
           f0=f1
15         f1=f
20      nm=k-1
        dk(0)=-sk(1)
        do 25 k=1,nm
25         dk(k)=-sk(k-1)-(k+1.0d0)/x*sk(k)
        return
        end
