        program mcjylv
c
c       ========================================================
c       purpose: this program computes bessel functions jv(z) 
c                and yv(z) and their derivatives with a large 
c                order and complex argument using subroutine
c                cjylv
c       input:   v --- order of jv(z) and yv(z)
c                z --- complex argument
c       output:  cbjv --- jv(z)
c                cdjv --- jv'(z)
c                cbyv --- yv(z)
c                cdyv --- yv'(z)
c       examples:
c                v = 100.00,    z = 4.00 + 2.00 i
c
c                jv(z) = -.6444792518-123 + .6619157435-123 i
c                jv'(z)= -.6251103777-122 + .1967638668-121 i
c                yv(z) =  .2403065353+121 + .2472039414+121 i
c                yv'(z)= -.7275814786+122 - .2533588851+122 i
c
c                v =100.5,     z = 4.00 + 2.00 i
c
c                jv(z) = -.1161315754-123 + .7390127781-124 i
c                jv'(z)= -.1588519437-122 + .2652227059-122 i
c                yv(z) =  .1941381412+122 + .1237578195+122 i
c                yv'(z)= -.5143285247+123 - .5320026773+122 i
c       ========================================================
c
        implicit double precision (v,x,y)
        implicit complex*16 (c,z)
        write(*,*)'please enter v,x and y ( z = x+iy )'
        read(*,*)v,x,y
        write(*,10)v,x,y
        z=cmplx(x,y)
        call cjylv(v,z,cbjv,cdjv,cbyv,cdyv)
        write(*,*)
        write(*,20)cbjv
        write(*,30)cdjv
        write(*,*)
        write(*,40)cbyv
        write(*,50)cdyv
10      format(8x,'v = ',f6.2,',    ','z =',f7.2,' + i ',f7.2)
20      format(8x,'jv(z) =',d17.10,' + i',d17.10)
30      format(8x,'jv''(z)=',d17.10,' + i',d17.10)
40      format(8x,'yv(z) =',d17.10,' + i',d17.10)
50      format(8x,'yv''(z)=',d17.10,' + i',d17.10)
        end


        subroutine cjylv(v,z,cbjv,cdjv,cbyv,cdyv)
c
c       ===================================================
c       purpose: compute bessel functions jv(z) and yv(z)
c                and their derivatives with a complex
c                argument and a large order
c       input:   v --- order of jv(z) and yv(z)
c                z --- complex argument
c       output:  cbjv --- jv(z)
c                cdjv --- jv'(z)
c                cbyv --- yv(z)
c                cdyv --- yv'(z)
c       routine called:
c                cjk to compute the expansion coefficients
c       ===================================================
c
        implicit double precision (a,b,d-h,o-y)
        implicit complex*16 (c,z)
        dimension cf(12),a(91)
        km=12
        call cjk(km,a)
        pi=3.141592653589793d0
        do 30 l=1,0,-1
           v0=v-l
           cws=cdsqrt(1.0d0-(z/v0)*(z/v0))
           ceta=cws+cdlog(z/v0/(1.0d0+cws))
           ct=1.0d0/cws
           ct2=ct*ct
           do 15 k=1,km
              l0=k*(k+1)/2+1
              lf=l0+k
              cf(k)=a(lf)
              do 10 i=lf-1,l0,-1
10               cf(k)=cf(k)*ct2+a(i)
15            cf(k)=cf(k)*ct**k
           vr=1.0d0/v0
           csj=(1.0d0,0.0d0)
           do 20 k=1,km
20            csj=csj+cf(k)*vr**k
           cbjv=cdsqrt(ct/(2.0d0*pi*v0))*cdexp(v0*ceta)*csj
           if (l.eq.1) cfj=cbjv
           csy=(1.0d0,0.0d0)
           do 25 k=1,km
25            csy=csy+(-1)**k*cf(k)*vr**k
           cbyv=-cdsqrt(2.0d0*ct/(pi*v0))*cdexp(-v0*ceta)*csy
           if (l.eq.1) cfy=cbyv
30      continue
        cdjv=-v/z*cbjv+cfj
        cdyv=-v/z*cbyv+cfy
        return
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

