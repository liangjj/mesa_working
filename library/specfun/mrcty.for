        program mrcty
c
c       =======================================================
c       purpose: this program computes the riccati-bessel 
c                functions of the second kind and their
c                derivatives using subroutine rcty
c       input:   x --- argument of riccati-bessel function
c                n --- order of yn(x)
c       output:  ry(n) --- xúyn(x)
c                dy(n) --- [xúyn(x)]'
c       example: x = 10.0
c                  n        xúyn(x)             [xúyn(x)]'
c                --------------------------------------------
c                  0     .8390715291d+00    -.5440211109d+00
c                  1     .6279282638d+00     .7762787027d+00
c                  2    -.6506930499d+00     .7580668738d+00
c                  3    -.9532747888d+00    -.3647106133d+00
c                  4    -.1659930220d-01    -.9466350679d+00
c                  5     .9383354168d+00    -.4857670106d+00
c       =======================================================
c
        implicit double precision (a-h,o-z)
        dimension ry(0:250),dy(0:250)
        write(*,*)'  please enter n and x '
        read(*,*)n,x
        write(*,30)n,x
        if (n.le.10) then
           ns=1
        else
           write(*,*)'  please enter order step ns'
           read(*,*)ns
        endif
        write(*,*)
        call rcty(n,x,nm,ry,dy)
        write(*,*)
        write(*,*)'  n        xúyn(x)             [xúyn(x)]'''
        write(*,*)'--------------------------------------------'
        do 10 k=0,nm,ns
           write(*,20)k,ry(k),dy(k)
10      continue
20      format(1x,i3,2d20.10)
30      format(3x,6hnmax =,i3,',    ',3hx =,f6.2)
        end


        subroutine rcty(n,x,nm,ry,dy)
c
c       ========================================================
c       purpose: compute riccati-bessel functions of the second
c                kind and their derivatives
c       input:   x --- argument of riccati-bessel function
c                n --- order of yn(x)
c       output:  ry(n) --- xúyn(x)
c                dy(n) --- [xúyn(x)]'
c                nm --- highest order computed
c       ========================================================
c
        implicit double precision (a-h,o-z)
        dimension ry(0:n),dy(0:n)
        nm=n
        if (x.lt.1.0d-60) then
           do 10 k=0,n
              ry(k)=-1.0d+300
10            dy(k)=1.0d+300
           ry(0)=-1.0d0
           dy(0)=0.0d0
           return
        endif
        ry(0)=-dcos(x)
        ry(1)=ry(0)/x-dsin(x)
        rf0=ry(0)
        rf1=ry(1)
        do 15 k=2,n
           rf2=(2.0d0*k-1.0d0)*rf1/x-rf0
           if (dabs(rf2).gt.1.0d+300) go to 20
           ry(k)=rf2
           rf0=rf1
15         rf1=rf2
20      nm=k-1
        dy(0)=dsin(x)
        do 25 k=1,nm
25         dy(k)=-k*ry(k)/x+ry(k-1)
        return
        end
