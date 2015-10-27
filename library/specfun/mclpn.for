        program mclpn
c
c       ==========================================================
c       purpose: this program computes the legendre polynomials 
c                pn(z) and pn'(z) for a complex argument using
c                subroutine clpn
c       input :  x --- real part of z
c                y --- imaginary part of z
c                n --- degree of pn(z), n = 0,1,...,n
c       output:  cpn(n) --- pn(z)
c                cpd(n) --- pn'(z)
c       example: z = 3.0 +2.0 i
c
c       n    re[pn(z)]     im[pn(z)]     re[pn'(z)]   im[pn'(z)]
c      -----------------------------------------------------------
c       0   .100000d+01   .000000d+00   .000000d+00   .000000d+00
c       1   .300000d+01   .200000d+01   .100000d+01   .000000d+00
c       2   .700000d+01   .180000d+02   .900000d+01   .600000d+01
c       3  -.270000d+02   .112000d+03   .360000d+02   .900000d+02
c       4  -.539000d+03   .480000d+03  -.180000d+03   .790000d+03
c       5  -.461700d+04   .562000d+03  -.481500d+04   .441000d+04
c       ==========================================================
c
        implicit double precision (x,y)
        implicit complex *16 (c,z)
        dimension cpn(0:100),cpd(0:100)
        write(*,*)'  please enter nmax, x and y (z=x+iy)'
        read(*,*)n,x,y
        write(*,30)x,y
        write(*,*)
        call clpn(n,x,y,cpn,cpd)
        write(*,*)'  n    re[pn(z)]     im[pn(z)]     re[pn''(z)]',
     &            '   im[pn''(z)]'
        write(*,*)' ---------------------------------------------',
     &            '--------------'
        do 10 k=0,n
10         write(*,20)k,cpn(k),cpd(k)
20      format(1x,i3,4d14.6)
30      format(3x,'x =',f5.1,',  ','y =',f5.1)
        end


        subroutine clpn(n,x,y,cpn,cpd)
c
c       ==================================================
c       purpose: compute legendre polynomials pn(z) and
c                their derivatives pn'(z) for a complex
c                argument
c       input :  x --- real part of z
c                y --- imaginary part of z
c                n --- degree of pn(z), n = 0,1,2,...
c       output:  cpn(n) --- pn(z)
c                cpd(n) --- pn'(z)
c       ==================================================
c
        implicit double precision (x,y)
        implicit complex *16 (c,z)
        dimension cpn(0:n),cpd(0:n)
        z=cmplx(x,y)
        cpn(0)=(1.0d0,0.0d0)
        cpn(1)=z
        cpd(0)=(0.0d0,0.0d0)
        cpd(1)=(1.0d0,0.0d0)
        cp0=(1.0d0,0.0d0)
        cp1=z
        do 10 k=2,n
           cpf=(2.0d0*k-1.0d0)/k*z*cp1-(k-1.0d0)/k*cp0
           cpn(k)=cpf
           if (dabs(x).eq.1.0d0.and.y.eq.0.0d0) then
              cpd(k)=0.5d0*x**(k+1)*k*(k+1.0d0)
           else
              cpd(k)=k*(cp1-z*cpf)/(1.0d0-z*z)
           endif
           cp0=cp1
10         cp1=cpf
        return
        end
