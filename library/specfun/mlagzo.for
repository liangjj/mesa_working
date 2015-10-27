        program mlagzo
c
c       ===========================================================
c       purpose : this program computes the zeros of laguerre
c                 polynomial ln(x) in the interval [0,ì] and the
c                 corresponding weighting coefficients for gauss-
c                 laguerre integration using subroutine lagzo
c       input :   n    --- order of the laguerre polynomial
c                 x(n) --- zeros of the laguerre polynomial
c                 w(n) --- corresponding weighting coefficients
c       ===========================================================
c
        implicit double precision (a-h,o-z)
        dimension x(100),w(100)
        write(*,*)'please enter the order of ln(x), n '
        read(*,*)n
        write(*,20)n
        call lagzo(n,x,w)
        write(*,*)'  nodes and weights for gauss-lagurre integration'
        write(*,*)
        write(*,*)'  i             xi                      wi'
        write(*,*)' -----------------------------------------',
     &            '------------'
        do 10 j=1,n
10         write(*,30)j,x(j),w(j)
20      format(1x,'n =',i3)
30      format(1x,i3,3x,d22.13,3x,d22.13)
        end


        subroutine lagzo(n,x,w)
c
c       =========================================================
c       purpose : compute the zeros of laguerre polynomial ln(x)
c                 in the interval [0,ì], and the corresponding
c                 weighting coefficients for gauss-laguerre
c                 integration
c       input :   n    --- order of the laguerre polynomial
c                 x(n) --- zeros of the laguerre polynomial
c                 w(n) --- corresponding weighting coefficients
c       =========================================================
c
        implicit double precision (a-h,o-z)
        dimension x(n),w(n)
        hn=1.0d0/n
        do 35 nr=1,n
           if (nr.eq.1) z=hn
           if (nr.gt.1) z=x(nr-1)+hn*nr**1.27
           it=0
10         it=it+1
           z0=z
           p=1.0d0
           do 15 i=1,nr-1
15            p=p*(z-x(i))
           f0=1.0d0
           f1=1.0d0-z
           do 20 k=2,n
              pf=((2.0d0*k-1.0d0-z)*f1-(k-1.0d0)*f0)/k
              pd=k/z*(pf-f1)
              f0=f1
20            f1=pf
           fd=pf/p
           q=0.0d0
           do 30 i=1,nr-1
              wp=1.0d0
              do 25 j=1,nr-1
                 if (j.eq.i) go to 25
                 wp=wp*(z-x(j))
25            continue
              q=q+wp
30         continue
           gd=(pd-q*fd)/p
           z=z-fd/gd
           if (it.le.40.and.dabs((z-z0)/z).gt.1.0d-15) go to 10
           x(nr)=z
           w(nr)=1.0d0/(z*pd*pd)
35      continue
        return
        end
