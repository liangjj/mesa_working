        program mpsi
c
c       ==================================================
c       purpose: this program computes the psi function
c                using subroutine psi
c       input :  x  --- argument of psi(x)
c       output:  ps --- psi(x)
c       examples:
c                   x          psi(x)
c                ------------------------
c                  .25      -4.227453533
c                  .50      -1.963510026
c                  .75      -1.085860880
c                 1.00       -.577215665
c                 1.25       -.227453533
c                 1.50        .036489974
c                 1.75        .247472454
c                 2.00        .422784335
c       ==================================================
c
        double precision x,ps
        write(*,*)'please enter x'
        read(*,*)x
        write(*,*)'    x          psi(x)'
        write(*,*)' ------------------------'
        call psi(x,ps)
        write(*,10)x,ps
10      format(1x,f6.2,f18.9)
        end


        subroutine psi(x,ps)
c
c       ======================================
c       purpose: compute the psi function
c       input :  x  --- argument of psi(x)
c       output:  ps --- psi(x)
c       ======================================
c
        implicit double precision (a-h,o-z)
        xa=dabs(x)
        pi=3.141592653589793d0
        el=.5772156649015329d0
        s=0.0d0
        if (x.eq.int(x).and.x.le.0.0) then
           ps=1.0d+300
           return
        else if (xa.eq.int(xa)) then
           n=xa
           do 10 k=1 ,n-1
10            s=s+1.0d0/k
           ps=-el+s
        else if (xa+.5.eq.int(xa+.5)) then
           n=xa-.5
           do 20 k=1,n
20            s=s+1.0/(2.0d0*k-1.0d0)
           ps=-el+2.0d0*s-1.386294361119891d0
        else
           if (xa.lt.10.0) then
              n=10-int(xa)
              do 30 k=0,n-1
30               s=s+1.0d0/(xa+k)
              xa=xa+n
           endif
           x2=1.0d0/(xa*xa)
           a1=-.8333333333333d-01
           a2=.83333333333333333d-02
           a3=-.39682539682539683d-02
           a4=.41666666666666667d-02
           a5=-.75757575757575758d-02
           a6=.21092796092796093d-01
           a7=-.83333333333333333d-01
           a8=.4432598039215686d0
           ps=dlog(xa)-.5d0/xa+x2*(((((((a8*x2+a7)*x2+
     &        a6)*x2+a5)*x2+a4)*x2+a3)*x2+a2)*x2+a1)
           ps=ps-s
        endif
        if (x.lt.0.0) ps=ps-pi*dcos(pi*x)/dsin(pi*x)-1.0d0/x
        return
        end
