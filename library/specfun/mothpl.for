        program mothpl
c
c       =========================================================
c       purpose: this program computes orthogonal polynomials: 
c                tn(x) or un(x) or ln(x) or hn(x), and their
c                derivatives using subroutine othpl
c       input :  kf --- function code
c                       kf=1 for chebyshev polynomial tn(x)
c                       kf=2 for chebyshev polynomial un(x)
c                       kf=3 for laguerre polynomial ln(x)
c                       kf=4 for hermite polynomial hn(x)
c                n ---  order of orthogonal polynomials
c                x ---  argument
c       output:  pl(n) --- tn(x) or un(x) or ln(x) or hn(x)
c                dpl(n)--- tn'(x) or un'(x) or ln'(x) or hn'(x)
c                          n = 0,1,2,...,n ( n ó 100 )
c       =========================================================
c
        implicit double precision (a-h,o-z)
        dimension pl(0:100),dpl(0:100)
        write(*,*)'kf,n,x = ?'
        read(*,*)kf,n,x
        write(*,20)kf,n,x
        write(*,*)
        call othpl(kf,n,x,pl,dpl)
        if (kf.eq.1) write(*,*)'  n          tn(x)            tn''(x)'
        if (kf.eq.2) write(*,*)'  n          un(x)            un''(x)'
        if (kf.eq.3) write(*,*)'  n          ln(x)            ln''(x)'
        if (kf.eq.4) write(*,*)'  n          hn(x)            hn''(x)'
        write(*,*)'-----------------------------------------'
        do 10 k=0,n
10         write(*,30)k,pl(k),dpl(k)
20      format(1x,3hkf=,i3,5x,5hnmax=,i3,5x,2hx=,f6.3)
30      format(1x,i3,2d18.8)
        end


        subroutine othpl(kf,n,x,pl,dpl)
c
c       ==========================================================
c       purpose: compute orthogonal polynomials: tn(x) or un(x),
c                or ln(x) or hn(x), and their derivatives
c       input :  kf --- function code
c                       kf=1 for chebyshev polynomial tn(x)
c                       kf=2 for chebyshev polynomial un(x)
c                       kf=3 for laguerre polynomial ln(x)
c                       kf=4 for hermite polynomial hn(x)
c                n ---  order of orthogonal polynomials
c                x ---  argument of orthogonal polynomials
c       output:  pl(n) --- tn(x) or un(x) or ln(x) or hn(x)
c                dpl(n)--- tn'(x) or un'(x) or ln'(x) or hn'(x)
c       =========================================================
c
        implicit double precision (a-h,o-z)
        dimension pl(0:n),dpl(0:n)
        a=2.0d0
        b=0.0d0
        c=1.0d0
        y0=1.0d0
        y1=2.0d0*x
        dy0=0.0d0
        dy1=2.0d0
        pl(0)=1.0d0
        pl(1)=2.0d0*x
        dpl(0)=0.0d0
        dpl(1)=2.0d0
        if (kf.eq.1) then
           y1=x
           dy1=1.0d0
           pl(1)=x
           dpl(1)=1.0d0
        else if (kf.eq.3) then
           y1=1.0d0-x
           dy1=-1.0d0
           pl(1)=1.0d0-x
           dpl(1)=-1.0d0
        endif
        do 10 k=2,n
           if (kf.eq.3) then
              a=-1.0d0/k
              b=2.0d0+a
              c=1.0d0+a
           else if (kf.eq.4) then
              c=2.0d0*(k-1.0d0)
           endif
           yn=(a*x+b)*y1-c*y0
           dyn=a*y1+(a*x+b)*dy1-c*dy0
           pl(k)=yn
           dpl(k)=dyn
           y0=y1
           y1=yn
           dy0=dy1
10         dy1=dyn
        return
        end
