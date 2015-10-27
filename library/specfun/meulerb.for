        program meulerb
c
c       ==========================================================
c       purpose: this program computes euler number en using
c                subroutine eulerb
c       example: compute euler number en for n = 0,2,...,10
c                computed results:
c
c                   n            en
c                 --------------------------
c                   0     .100000000000d+01
c                   2    -.100000000000d+01
c                   4     .500000000000d+01
c                   6    -.610000000000d+02
c                   8     .138500000000d+04
c                  10    -.505210000000d+05
c       ==========================================================
c
        double precision e
        dimension e(0:200)
        write(*,*)'  please enter nmax '
        read(*,*)n
        call eulerb(n,e)
        write(*,*)'   n            en'
        write(*,*)' --------------------------'
        do 10 k=0,n,2
10         write(*,20)k,e(k)
20      format(2x,i3,d22.12)
        end


        subroutine eulerb(n,en)
c
c       ======================================
c       purpose: compute euler number en
c       input :  n --- serial number
c       output:  en(n) --- en
c       ======================================
c
        implicit double precision (a-h,o-z)
        dimension en(0:n)
        hpi=2.0d0/3.141592653589793d0
        en(0)=1.0d0
        en(2)=-1.0d0
        r1=-4.0d0*hpi**3
        do 20 m=4,n,2
           r1=-r1*(m-1)*m*hpi*hpi
           r2=1.0d0
           isgn=1.0d0
           do 10 k=3,1000,2
              isgn=-isgn
              s=(1.0d0/k)**(m+1)
              r2=r2+isgn*s
              if (s.lt.1.0d-15) goto 20
10         continue
20         en(m)=r1*r2
        return
        end
