        program mciklv
c
c       =========================================================
c       purpose: this program computes modified bessel functions 
c                iv(z) and kv(z) and their derivatives for a  
c                large order and a complex argument using
c                subroutine ciklv
c       input:   v --- order of iv(z) and kv(z)
c                z --- complex argument
c       output:  cbiv --- iv(z)
c                cdiv --- iv'(z)
c                cbkv --- kv(z)
c                cdkv --- kv'(z)
c       examples:
c                v =100.00,    z =   4.00 + i   2.00
c
c       iv(z) = -.7373606617-123 + .6461109082-123 i
c       iv'(z)= -.8307094243-122 + .2030132500-121 i
c       kv(z) = -.3836166007+121 - .3356017795+121 i
c       kv'(z)=  .1103271276+123 + .2886519240+122 i
c
c                v =100.50,    z =   4.00 + i   2.00
c       iv(z) = -.1289940051-123 + .6845756182-124 i
c       iv'(z)= -.1907996261-122 + .2672465997-122 i
c       kv(z) = -.3008779281+122 - .1593719779+122 i
c       kv'(z)=  .7653781978+123 + .1857772148+122 i
c       =========================================================
c
        implicit double precision (v,x,y)
        implicit complex*16 (c,z)
        write(*,*)'please enter v,x,y ( z = x+iy )'
        read(*,*)v,x,y
        write(*,10)v,x,y
        z=cmplx(x,y)
        call ciklv(v,z,cbiv,cdiv,cbkv,cdkv)
        write(*,*)
        write(*,20)cbiv
        write(*,30)cdiv
        write(*,*)
        write(*,40)cbkv
        write(*,50)cdkv
10      format(8x,'v =',f6.2,',    ','z =',f7.2,' + i',f7.2)
20      format(8x,'iv(z) =',d17.10,' + i ',d17.10)
30      format(8x,'iv''(z)=',d17.10,' + i ',d17.10)
40      format(8x,'kv(z) =',d17.10,' + i ',d17.10)
50      format(8x,'kv''(z)=',d17.10,' + i ',d17.10)
        end



        subroutine ciklv(v,z,cbiv,cdiv,cbkv,cdkv)
c
c       =====================================================
c       purpose: compute modified bessel functions iv(z) and
c                kv(z) and their derivatives with a complex
c                argument and a large order
c       input:   v --- order of iv(z) and kv(z)
c                z --- complex argument
c       output:  cbiv --- iv(z)
c                cdiv --- iv'(z)
c                cbkv --- kv(z)
c                cdkv --- kv'(z)
c       routine called:
c                cjk to compute the expansion coefficients
c       ====================================================
c
        implicit double precision (a,b,d-h,o-y)
        implicit complex*16 (c,z)
        dimension cf(12),a(91)
        pi=3.141592653589793d0
        km=12
        call cjk(km,a)
        do 30 l=1,0,-1
           v0=v-l
           cws=cdsqrt(1.0d0+(z/v0)*(z/v0))
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
           csi=(1.0d0,0.0d0)
           do 20 k=1,km
20            csi=csi+cf(k)*vr**k
           cbiv=cdsqrt(ct/(2.0d0*pi*v0))*cdexp(v0*ceta)*csi
           if (l.eq.1) cfi=cbiv
           csk=(1.0d0,0.0d0)
           do 25 k=1,km
25            csk=csk+(-1)**k*cf(k)*vr**k
           cbkv=cdsqrt(pi*ct/(2.0d0*v0))*cdexp(-v0*ceta)*csk
           if (l.eq.1) cfk=cbkv
30      continue
        cdiv=cfi-v/z*cbiv
        cdkv=-cfk-v/z*cbkv
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

