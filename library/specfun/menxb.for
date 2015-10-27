        program menxb
c
c       =========================================================
c       purpose: this program computes the exponential integral 
c                en(x) using subroutine enxb
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
        call enxb(n,x,en)
        do 10 k=0,n
           write(*,30)k,en(k)
10      continue
20      format(5x,i3,',   ','x=',f5.1)
30      format(2x,i3,d18.8)
        end


        subroutine enxb(n,x,en)
c
c       ===============================================
c       purpose: compute exponential integral en(x)
c       input :  x --- argument of en(x)
c                n --- order of en(x)  (n = 0,1,2,...)
c       output:  en(n) --- en(x)
c       ===============================================
c
        implicit double precision (a-h,o-z)
        dimension en(0:n)
        if (x.eq.0.0) then
           en(0)=1.0d+300
           en(1)=1.0d+300
           do 10 k=2,n
10            en(k)=1.0d0/(k-1.0)
           return
        else if (x.le.1.0) then
           en(0)=dexp(-x)/x
           do 40 l=1,n
              rp=1.0d0
              do 15 j=1,l-1
15               rp=-rp*x/j
              ps=-0.5772156649015328d0
              do 20 m=1,l-1
20               ps=ps+1.0d0/m
              ens=rp*(-dlog(x)+ps)
              s=0.0d0
              do 30 m=0,20
                 if (m.eq.l-1) go to 30
                 r=1.0d0
                 do 25 j=1,m
25                  r=-r*x/j
                 s=s+r/(m-l+1.0d0)
                 if (dabs(s-s0).lt.dabs(s)*1.0d-15) go to 35
                 s0=s
30            continue
35            en(l)=ens-s
40         continue
        else
           en(0)=dexp(-x)/x
           m=15+int(100.0/x)
           do 50 l=1,n
              t0=0.0d0
              do 45 k=m,1,-1
45               t0=(l+k-1.0d0)/(1.0d0+k/(x+t0))
              t=1.0d0/(x+t0)
50            en(l)=dexp(-x)*t
        endif
        end
