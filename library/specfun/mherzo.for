        program mherzo
c
c       ===========================================================
c       purpose : this program computes the zeros of hermite 
c                 polynomial ln(x) in the interval [-ì,ì] and the
c                 corresponding weighting coefficients for gauss-
c                 hermite integration using subroutine herzo
c       input :   n    --- order of the hermite polynomial
c                 x(n) --- zeros of the hermite polynomial
c                 w(n) --- corresponding weighting coefficients
c       ===========================================================
c
        implicit double precision (a-h,o-z)
        dimension x(100),w(100)
        write(*,*)'please enter the order of hn(x), n '
        read(*,*)n
        write(*,20)n
        call herzo(n,x,w)
        write(*,*)'  nodes and weights for gauss-hermite integration'
        write(*,*)
        write(*,*)'  i             xi                      wi'
        write(*,*)' -----------------------------------------',
     &            '------------'
        do 10 j=1,n
10         write(*,30)j,x(j),w(j)
20      format(1x,'n =',i3)
30      format(1x,i3,3x,d22.13,3x,d22.13)
        end


        subroutine herzo(n,x,w)
c
c       ========================================================
c       purpose : compute the zeros of hermite polynomial ln(x)
c                 in the interval [-ì,ì], and the corresponding
c                 weighting coefficients for gauss-hermite
c                 integration
c       input :   n    --- order of the hermite polynomial
c                 x(n) --- zeros of the hermite polynomial
c                 w(n) --- corresponding weighting coefficients
c       ========================================================
c
        implicit double precision (a-h,o-z)
        dimension x(n),w(n)
        hn=1.0d0/n
        zl=-1.1611d0+1.46d0*n**0.5
        do 40 nr=1,n/2
           if (nr.eq.1) z=zl
           if (nr.ne.1) z=z-hn*(n/2+1-nr)
           it=0
10         it=it+1
           z0=z
           f0=1.0d0
           f1=2.0d0*z
           do 15 k=2,n
              hf=2.0d0*z*f1-2.0d0*(k-1.0d0)*f0
              hd=2.0d0*k*f1
              f0=f1
15            f1=hf
           p=1.0d0
           do 20 i=1,nr-1
20            p=p*(z-x(i))
           fd=hf/p
           q=0.0d0
           do 30 i=1,nr-1
              wp=1.0d0
              do 25 j=1,nr-1
                 if (j.eq.i) go to 25
                 wp=wp*(z-x(j))
25            continue
30            q=q+wp
           gd=(hd-q*fd)/p
           z=z-fd/gd
           if (it.le.40.and.dabs((z-z0)/z).gt.1.0d-15) go to 10
           x(nr)=z
           x(n+1-nr)=-z
           r=1.0d0
           do 35 k=1,n
35            r=2.0d0*r*k
           w(nr)=3.544907701811d0*r/(hd*hd)
40         w(n+1-nr)=w(nr)
        if (n.ne.2*int(n/2)) then
           r1=1.0d0
           r2=1.0d0
           do 45 j=1,n
              r1=2.0d0*r1*j
              if (j.ge.(n+1)/2) r2=r2*j
45         continue
           w(n/2+1)=0.88622692545276d0*r1/(r2*r2)
           x(n/2+1)=0.0d0
        endif
        return
        end
