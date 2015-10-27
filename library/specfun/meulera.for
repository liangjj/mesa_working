        program meulera
c
c       ==========================================================
c       purpose: this program computes euler number en using
c                subroutine eulera
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
        call eulera(n,e)
        write(*,*)'   n            en'
        write(*,*)' --------------------------'
        do 10 k=0,n,2
10         write(*,20)k,e(k)
20      format(2x,i3,d22.12)
        end


        subroutine eulera(n,en)
c
c       ======================================
c       purpose: compute euler number en
c       input :  n --- serial number
c       output:  en(n) --- en
c       ======================================
c
        implicit double precision (a-h,o-z)
        dimension en(0:n)
        en(0)=1.0d0
        do 30 m=1,n/2
           s=1.0d0
           do 20 k=1,m-1
              r=1.0d0
              do 10 j=1,2*k
10               r=r*(2.0d0*m-2.0d0*k+j)/j
20            s=s+r*en(2*k)
30         en(2*m)=-s
        return
        end
