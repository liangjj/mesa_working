        program mcerzo
c
c       ===============================================================
c       purpose : this program evaluates the complex zeros of error 
c                 function erf(z) using subroutine cerzo
c       input:    nt --- total number of zeros
c       example:  nt = 10
c
c    n     complex zeros of erf(z)     n     complex zeros of erf(z)
c   -------------------------------------------------------------------
c    1   1.450616163 + i 1.880943000   6   4.158998400 + i 4.435571444
c    2   2.244659274 + i 2.616575141   7   4.516319400 + i 4.780447644
c    3   2.839741047 + i 3.175628100   8   4.847970309 + i 5.101588043
c    4   3.335460735 + i 3.646174376   9   5.158767908 + i 5.403332643
c    5   3.769005567 + i 4.060697234  10   5.452192201 + i 5.688837437
c       ===============================================================
c
        implicit double precision (e,p,w)
        implicit complex *16 (c,z)
        dimension zo(100)
        write(*,*)'please enter nt '
        read(*,*)nt
        write(*,20)nt
        call cerzo(nt,zo)
        write(*,*)'  *****    please wait !    *****'
        write(*,*)
        write(*,*)'  n        complex zeros of erf(z)'
        write(*,*)'-------------------------------------'
        do 10 i=1,nt
10         write(*,30) i,zo(i)
20      format(2x,'nt=',i3)
30      format(1x,i3,2x,f13.8,2x,2h+i,f13.8)
        end


        subroutine cerzo(nt,zo)
c
c       ===============================================================
c       purpose : evaluate the complex zeros of error function erf(z)
c                 using the modified newton's iteration method
c       input :   nt --- total number of zeros
c       output:   zo(l) --- l-th zero of erf(z), l=1,2,...,nt
c       routine called: cerf for computing erf(z) and erf'(z)
c       ===============================================================
c
        implicit double precision (e,p,w)
        implicit complex *16 (c,z)
        dimension zo(nt)
        pi=3.141592653589793d0
        do 35 nr=1,nt
           pu=dsqrt(pi*(4.0d0*nr-0.5d0))
           pv=pi*dsqrt(2.0d0*nr-0.25d0)
           px=0.5*pu-0.5*dlog(pv)/pu
           py=0.5*pu+0.5*dlog(pv)/pu
           z=cmplx(px,py)
           it=0
15         it=it+1
           call cerf(z,zf,zd)
           zp=(1.0d0,0.0d0)
           do 20 i=1,nr-1
20            zp=zp*(z-zo(i))
           zfd=zf/zp
           zq=(0.0d0,0.0d0)
           do 30 i=1,nr-1
              zw=(1.0d0,0.0d0)
              do 25 j=1,nr-1
                 if (j.eq.i) go to 25
                 zw=zw*(z-zo(j))
25            continue
30            zq=zq+zw
           zgd=(zd-zq*zfd)/zp
           z=z-zfd/zgd
           w0=w
           w=cdabs(z)
           if (it.le.50.and.dabs((w-w0)/w).gt.1.0d-11) go to 15
35         zo(nr)=z
        return
        end


        subroutine cerf(z,cer,cder)
c
c       ==========================================================
c       purpose: compute complex error function erf(z) & erf'(z)
c       input:   z   --- complex argument of erf(z)
c                x   --- real part of z
c                y   --- imaginary part of z
c       output:  cer --- erf(z)
c                cder --- erf'(z)
c       ==========================================================
        implicit double precision (a-h,o-z)
        complex *16 z,cer,cder
        eps=1.0d-12
        pi=3.141592653589793d0
        x=real(z)
        y=dimag(z)
        x2=x*x
        if (x.le.3.5d0) then
           er=1.0d0
           r=1.0d0
           do 10 k=1,100
              r=r*x2/(k+0.5d0)
              er=er+r
              if (dabs(er-w).le.eps*dabs(er)) go to 15
10            w=er
15         c0=2.0d0/dsqrt(pi)*x*dexp(-x2)
           er0=c0*er
        else
           er=1.0d0
           r=1.0d0
           do 20 k=1,12
              r=-r*(k-0.5d0)/x2
20            er=er+r
           c0=dexp(-x2)/(x*dsqrt(pi))
           er0=1.0d0-c0*er
        endif
        if (y.eq.0.0d0) then
           err=er0
           eri=0.0d0
        else
           cs=dcos(2.0d0*x*y)
           ss=dsin(2.0d0*x*y)
           er1=dexp(-x2)*(1.0d0-cs)/(2.0d0*pi*x)
           ei1=dexp(-x2)*ss/(2.0d0*pi*x)
           er2=0.0d0
           do 25 n=1,100
              er2=er2+dexp(-.25d0*n*n)/(n*n+4.0d0*x2)*(2.0d0*x
     &            -2.0d0*x*dcosh(n*y)*cs+n*dsinh(n*y)*ss)
              if (dabs((er2-w1)/er2).lt.eps) go to 30
25            w1=er2
30         c0=2.0d0*dexp(-x2)/pi
           err=er0+er1+c0*er2
           ei2=0.0d0
           do 35 n=1,100
              ei2=ei2+dexp(-.25d0*n*n)/(n*n+4.0d0*x2)*(2.0d0*x
     &            *dcosh(n*y)*ss+n*dsinh(n*y)*cs)
              if (dabs((ei2-w2)/ei2).lt.eps) go to 40
35            w2=ei2
40         eri=ei1+c0*ei2
        endif
        cer=cmplx(err,eri)
        cder=2.0d0/dsqrt(pi)*cdexp(-z*z)
        return
        end
