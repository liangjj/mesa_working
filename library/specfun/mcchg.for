        program mcchg
c
c       ===========================================================
c       purpose: this program computes confluent hypergeometric 
c                function m(a,b,z) with real parameters a, b, and 
c                a complex argument z using subroutine cchg
c       input :  a --- parameter
c                b --- parameter
c                z --- complex argument
c       output:  chg --- m(a,b,z)
c       examples:
c          a      b        z        re[m(a,b,z)]   im[m(a,b,z)]
c         -------------------------------------------------------
c         3.3   4.25    10 + 0i    .61677489d+04    0
c         3.3   4.25    25 + 0i    .95781835d+10  -.15738228d-03
c         3.3   4.25     3 -  i    .75828716d+01  -.86815474d+01
c         3.3   4.25    15 +10i   -.58313765d+06  -.48195426d+05
c       ===========================================================
c
        implicit double precision (a,b,x,y)
        implicit complex *16 (c,z)
        write(*,*)'please enter a, b, x and y (z=x+iy) '
        read(*,*)a,b,x,y
        write(*,20)a,b,x,y
        z=cmplx(x,y)
        call cchg(a,b,z,chg)
        write(*,10)chg
10      format(10x,'m(a,b,z) =',d18.8,' + i ',d18.8)
20      format(1x,'a =',f5.1,',  ','b =',f5.1,',  ','x =',f5.1,
     &         ',  ','y =',f5.1)
        end


        subroutine cchg(a,b,z,chg)
c
c       ===================================================
c       purpose: compute confluent hypergeometric function
c                m(a,b,z) with real parameters a, b and a
c                complex argument z
c       input :  a --- parameter
c                b --- parameter
c                z --- complex argument
c       output:  chg --- m(a,b,z)
c       routine called: gamma for computing gamma function
c       ===================================================
c
        implicit double precision (a,b,d-h,o-y)
        implicit complex *16 (c,z)
        pi=3.141592653589793d0
        ci=(0.0d0,1.0d0)
        a0=a
        a1=a
        z0=z
        if (b.eq.0.0.or.b.eq.-int(abs(b))) then
           chg=(1.0d+300,0.0d0)
        else if (a.eq.0.0d0.or.z.eq.0.0d0) then
           chg=(1.0d0,0.0d0)
        else if (a.eq.-1.0d0) then
           chg=1.0d0-z/b
        else if (a.eq.b) then
           chg=cdexp(z)
        else if (a-b.eq.1.0d0) then
           chg=(1.0d0+z/b)*cdexp(z)
        else if (a.eq.1.0d0.and.b.eq.2.0d0) then
           chg=(cdexp(z)-1.0d0)/z
        else if (a.eq.int(a).and.a.lt.0.0d0) then
           m=int(-a)
           cr=(1.0d0,0.0d0)
           chg=(1.0d0,0.0d0)
           do 10 k=1,m
              cr=cr*(a+k-1.0d0)/k/(b+k-1.0d0)*z
10            chg=chg+cr
        else
           x0=real(z)
           if (x0.lt.0.0d0) then
              a=b-a
              a0=a
              z=-z
           endif
           if (a.lt.2.0d0) nl=0
           if (a.ge.2.0d0) then
              nl=1
              la=int(a)
              a=a-la-1.0d0
           endif
           do 30 n=0,nl
              if (a0.ge.2.0d0) a=a+1.0d0
              if (cdabs(z).lt.20.0d0+abs(b).or.a.lt.0.0d0) then
                 chg=(1.0d0,0.0d0)
                 crg=(1.0d0,0.0d0)
                 do 15 j=1,500
                    crg=crg*(a+j-1.0d0)/(j*(b+j-1.0d0))*z
                    chg=chg+crg
                    if (cdabs((chg-chw)/chg).lt.1.d-15) go to 25
                    chw=chg
15               continue
              else
                 call gamma(a,g1)
                 call gamma(b,g2)
                 ba=b-a
                 call gamma(ba,g3)
                 cs1=(1.0d0,0.0d0)
                 cs2=(1.0d0,0.0d0)
                 cr1=(1.0d0,0.0d0)
                 cr2=(1.0d0,0.0d0)
                 do 20 i=1,8
                    cr1=-cr1*(a+i-1.0d0)*(a-b+i)/(z*i)
                    cr2=cr2*(b-a+i-1.0d0)*(i-a)/(z*i)
                    cs1=cs1+cr1
20                  cs2=cs2+cr2
                 x=real(z)
                 y=dimag(z)
                 if (x.eq.0.0.and.y.ge.0.0) then
                    phi=0.5d0*pi
                 else if (x.eq.0.0.and.y.le.0.0) then
                    phi=-0.5d0*pi
                 else
                    phi=datan(y/x)
                 endif
                 if (phi.gt.-0.5*pi.and.phi.lt.1.5*pi) ns=1
                 if (phi.gt.-1.5*pi.and.phi.le.-0.5*pi) ns=-1
                 cfac=cdexp(ns*ci*pi*a)
                 if (y.eq.0.0d0) cfac=dcos(pi*a)
                 chg1=g2/g3*z**(-a)*cfac*cs1
                 chg2=g2/g1*cdexp(z)*z**(a-b)*cs2
                 chg=chg1+chg2
              endif
25            if (n.eq.0) cy0=chg
              if (n.eq.1) cy1=chg
30         continue
           if (a0.ge.2.0d0) then
              do 35 i=1,la-1
                 chg=((2.0d0*a-b+z)*cy1+(b-a)*cy0)/a
                 cy0=cy1
                 cy1=chg
35               a=a+1.0d0
           endif
           if (x0.lt.0.0d0) chg=chg*cdexp(-z)
        endif
        a=a1
        z=z0
        return
        end


        subroutine gamma(x,ga)
c
c       ==================================================
c       purpose: compute gamma function ג(x)
c       input :  x  --- argument of ג(x)
c                       ( x is not equal to 0,-1,-2,תתת)
c       output:  ga --- ג(x)
c       ==================================================
c
        implicit double precision (a-h,o-z)
        dimension g(26)
        pi=3.141592653589793d0
        if (x.eq.int(x)) then
           if (x.gt.0.0d0) then
              ga=1.0d0
              m1=x-1
              do 10 k=2,m1
10               ga=ga*k
           else
              ga=1.0d+300
           endif
        else
           if (dabs(x).gt.1.0d0) then
              z=dabs(x)
              m=int(z)
              r=1.0d0
              do 15 k=1,m
15               r=r*(z-k)
              z=z-m
           else
              z=x
           endif
           data g/1.0d0,0.5772156649015329d0,
     &          -0.6558780715202538d0, -0.420026350340952d-1,
     &          0.1665386113822915d0,-.421977345555443d-1,
     &          -.96219715278770d-2, .72189432466630d-2,
     &          -.11651675918591d-2, -.2152416741149d-3,
     &          .1280502823882d-3, -.201348547807d-4,
     &          -.12504934821d-5, .11330272320d-5,
     &          -.2056338417d-6, .61160950d-8,
     &          .50020075d-8, -.11812746d-8,
     &          .1043427d-9, .77823d-11,
     &          -.36968d-11, .51d-12,
     &          -.206d-13, -.54d-14, .14d-14, .1d-15/
           gr=g(26)
           do 20 k=25,1,-1
20            gr=gr*z+g(k)
           ga=1.0d0/(gr*z)
           if (dabs(x).gt.1.0d0) then
              ga=ga*r
              if (x.lt.0.0d0) ga=-pi/(x*ga*dsin(pi*x))
           endif
        endif
        return
        end
