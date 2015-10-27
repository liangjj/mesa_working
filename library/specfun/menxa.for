        program menxa
c
c       =========================================================
c       purpose: this program computes the exponential integral 
c                en(x) using subroutine enxa
c       example: x = 10.0
c                   n         en(x)
c                 ----------------------
c                   0     .45399930d-05
c                   1     .41569689d-05
c                   2     .38302405d-05
c                   3     .35487626d-05
c                   4     .33041014d-05
c                   5     .30897289d-05
c       =========================================================
c
        implicit double precision (a-h,o-z)
        dimension en(0:100)
        write(*,*)'please enter n and x '
        read(*,*)n,x
        write(*,20)n,x
        write(*,*)
        write(*,*)'   n         en(x)'
        write(*,*)' ----------------------'
        call enxa(n,x,en)
        do 10 k=0,n
           write(*,30)k,en(k)
10      continue
20      format(5x,i3,',   ','x=',f5.1)
30      format(2x,i3,d18.8)
        end


        subroutine enxa(n,x,en)
c
c       ============================================
c       purpose: compute exponential integral en(x)
c       input :  x --- argument of en(x) ( x ó 20 )
c                n --- order of en(x)
c       output:  en(n) --- en(x)
c       routine called: e1xb for computing e1(x)
c       ============================================
c
        implicit double precision (a-h,o-z)
        dimension en(0:n)
        en(0)=dexp(-x)/x
        call e1xb(x,e1)
        en(1)=e1
        do 10 k=2,n
           ek=(dexp(-x)-x*e1)/(k-1.0d0)
           en(k)=ek
10         e1=ek
        return
        end


        subroutine e1xb(x,e1)
c
c       ============================================
c       purpose: compute exponential integral e1(x)
c       input :  x  --- argument of e1(x)
c       output:  e1 --- e1(x)
c       ============================================
c
        implicit double precision (a-h,o-z)
        if (x.eq.0.0) then
           e1=1.0d+300
        else if (x.le.1.0) then
           e1=1.0d0
           r=1.0d0
           do 10 k=1,25
              r=-r*k*x/(k+1.0d0)**2
              e1=e1+r
              if (dabs(r).le.dabs(e1)*1.0d-15) go to 15
10         continue
15         ga=0.5772156649015328d0
           e1=-ga-dlog(x)+x*e1
        else
           m=20+int(80.0/x)
           t0=0.0d0
           do 20 k=m,1,-1
              t0=k/(1.0d0+k/(x+t0))
20         continue
           t=1.0d0/(x+t0)
           e1=dexp(-x)*t
        endif
        return
        end
