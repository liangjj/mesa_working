        program mbernoa
c
c       ===========================================================
c       purpose: this program computes bernoulli number bn using 
c                subroutine bernoa
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
        call bernoa(n,b)
        write(*,*)'   n            bn'
        write(*,*)' --------------------------'
        write(*,20)0,b(0)
        write(*,20)1,b(1)
        do 10 k=2,n,2
10         write(*,20)k,b(k)
20      format(2x,i3,d22.12)
        end


        subroutine bernoa(n,bn)
c
c       ======================================
c       purpose: compute bernoulli number bn
c       input :  n --- serial number
c       output:  bn(n) --- bn
c       ======================================
c
        implicit double precision (a-h,o-z)
        dimension bn(0:n)
        bn(0)=1.0d0
        bn(1)=-0.5d0
        do 30 m=2,n
           s=-(1.0d0/(m+1.0d0)-0.5d0)
           do 20 k=2,m-1
              r=1.0d0
              do 10 j=2,k
10               r=r*(j+m-k)/j
20            s=s-r*bn(k)
30         bn(m)=s
        do 40 m=3,n,2
40         bn(m)=0.0d0
        return
        end
