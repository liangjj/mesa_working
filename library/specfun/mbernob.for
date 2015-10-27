        program mbernob
c
c       ===========================================================
c       purpose: this program computes bernoulli number bn using 
c                subroutine bernob
c       example: compute bernouli number bn for n = 0,1,...,10
c                computed results:
c
c                   n            bn
c                 --------------------------
c                   0     .100000000000d+01
c                   1    -.500000000000d+00
c                   2     .166666666667d+00
c                   4    -.333333333333d-01
c                   6     .238095238095d-01
c                   8    -.333333333333d-01
c                  10     .757575757576d-01
c       ===========================================================
c
        double precision b
        dimension b(0:200)
        write(*,*)'  please enter nmax'
        read(*,*)n
        call bernob(n,b)
        write(*,*)'   n            bn'
        write(*,*)' --------------------------'
        write(*,20)0,b(0)
        write(*,20)1,b(1)
        do 10 k=2,n,2
10         write(*,20)k,b(k)
20      format(2x,i3,d22.12)
        end


        subroutine bernob(n,bn)
c
c       ======================================
c       purpose: compute bernoulli number bn
c       input :  n --- serial number
c       output:  bn(n) --- bn
c       ======================================
c
        implicit double precision (a-h,o-z)
        dimension bn(0:n)
        tpi=6.283185307179586d0 
        bn(0)=1.0d0
        bn(1)=-0.5d0
        bn(2)=1.0d0/6.0d0
        r1=(2.0d0/tpi)**2
        do 20 m=4,n,2
           r1=-r1*(m-1)*m/(tpi*tpi)
           r2=1.0d0
           do 10 k=2,10000
              s=(1.0d0/k)**m
              r2=r2+s
              if (s.lt.1.0d-15) goto 20
10         continue
20         bn(m)=r1*r2
        return
        end
