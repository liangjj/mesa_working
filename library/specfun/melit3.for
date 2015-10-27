        program melit3
c
c       ==========================================================
c       purpose: this program computes the elliptic integral of 
c                the third kind using subroutine elit3
c       input :  phi --- argument ( in degrees )
c                 k  --- modulus   ( 0 ó k ó 1 )
c                 c  --- parameter ( 0 ó c ó 1 )
c       output:  el3 ÄÄÄ value of the elliptic integral of the
c                        third kind
c       ==========================================================
c
        implicit double precision (a-h,o-z)
        write(*,*)'please enter phi, k and c '
        read(*,*)phi,hk,c
        call elit3(phi,hk,c,el3)
        write(*,10)el3
10      format(1x,'el3=',f12.8)
        end


        subroutine elit3(phi,hk,c,el3)
c
c       =========================================================
c       purpose: compute the elliptic integral of the third kind
c                using gauss-legendre quadrature
c       input :  phi --- argument ( in degrees )
c                 k  --- modulus   ( 0 ó k ó 1.0 )
c                 c  --- parameter ( 0 ó c ó 1.0 )
c       output:  el3 --- value of the elliptic integral of the
c                        third kind
c       =========================================================
c
        implicit double precision (a-h,o-z)
        dimension t(10),w(10)
        logical lb1,lb2
        data t/.9931285991850949,.9639719272779138,
     &         .9122344282513259,.8391169718222188,
     &         .7463319064601508,.6360536807265150,
     &         .5108670019508271,.3737060887154195,
     &         .2277858511416451,.7652652113349734d-1/
        data w/.1761400713915212d-1,.4060142980038694d-1,
     &         .6267204833410907d-1,.8327674157670475d-1,
     &         .1019301198172404,.1181945319615184,
     &         .1316886384491766,.1420961093183820,
     &         .1491729864726037,.1527533871307258/
        lb1=hk.eq.1.0d0.and.dabs(phi-90.0).le.1.0d-8
        lb2=c.eq.1.0d0.and.dabs(phi-90.0).le.1.0d-8
        if (lb1.or.lb2) then
            el3=1.0d+300
            return
        endif
        c1=0.87266462599716d-2*phi
        c2=c1
        el3=0.0d0
        do 10 i=1,10
           c0=c2*t(i)
           t1=c1+c0
           t2=c1-c0
           f1=1.0d0/((1.0d0-c*dsin(t1)*dsin(t1))*
     &              dsqrt(1.0d0-hk*hk*dsin(t1)*dsin(t1)))
           f2=1.0d0/((1.0d0-c*dsin(t2)*dsin(t2))*
     &              dsqrt(1.0d0-hk*hk*dsin(t2)*dsin(t2)))
10         el3=el3+w(i)*(f1+f2)
        el3=c1*el3
        return
        end
