        program mcjk
c
c       ============================================================
c       purpose: this program computes the expansion coefficients  
c                for the asymptotic expansion of bessel functions  
c                with large orders using subroutine cjk
c       input :  km   --- maximum k
c       output:  a(l) --- cj(k) where j and k are related to l by
c                         l=j+1+[k*(k+1)]/2; j,k=0,1,2,...,km
c       ============================================================
c
        implicit double precision (a-h,o-z)
        dimension a(231)
        write(*,*)'please enter km ( ó 20 )'
        read(*,*)km
        lm=km+1+(km*(km+1))/2
        call cjk(km,a)
        do 10 k=1,lm
10         write(*,15)k,a(k)
15      format(1x,i3,d25.14)
        end


        subroutine cjk(km,a)
c
c       ========================================================
c       purpose: compute the expansion coefficients for the
c                asymptotic expansion of bessel functions 
c                with large orders
c       input :  km   --- maximum k
c       output:  a(l) --- cj(k) where j and k are related to l 
c                         by l=j+1+[k*(k+1)]/2; j,k=0,1,...,km
c       ========================================================
c
        implicit double precision (a-h,o-z)
        dimension a(*)
        a(1)=1.0d0
        f0=1.0d0
        g0=1.0d0
        do 10 k=0,km-1
           l1=(k+1)*(k+2)/2+1
           l2=(k+1)*(k+2)/2+k+2
           f=(0.5d0*k+0.125d0/(k+1))*f0
           g=-(1.5d0*k+0.625d0/(3.0*(k+1.0d0)))*g0
           a(l1)=f
           a(l2)=g
           f0=f
10         g0=g
        do 15 k=1,km-1
           do 15 j=1,k
              l3=k*(k+1)/2+j+1
              l4=(k+1)*(k+2)/2+j+1
              a(l4)=(j+0.5d0*k+0.125d0/(2.0*j+k+1.0))*a(l3)
     &             -(j+0.5d0*k-1.0+0.625d0/(2.0*j+k+1.0))*a(l3-1)
15         continue
        return
        end
