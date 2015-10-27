        program me1z
c
c       =========================================================
c       purpose: this program computes the complex exponential 
c                integral e1(z) using subroutine e1z
c       example:
c                     z            re[e1(z)]       im[e1(z)]
c                -----------------------------------------------
c                 3.0    2.0    -.90959209d-02   -.69001793d-02
c                 3.0   -2.0    -.90959209d-02    .69001793d-02
c                -3.0    2.0    -.28074891d+01    .59603353d+01
c                -3.0   -2.0    -.28074891d+01   -.59603353d+01
c                25.0   10.0    -.29302080d-12    .40391222d-12
c                25.0  -10.0    -.29302080d-12   -.40391222d-12
c               -25.0   10.0     .27279957d+10   -.49430610d+09
c               -25.0  -10.0     .27279957d+10    .49430610d+09
c       =========================================================
c
        implicit complex*16 (c,z)
        implicit double precision (d-h,o-y)
        write(*,*)'please enter x and y ( z =x+iy ) '
        read(*,*)x,y
        z=cmplx(x,y)
        call e1z(z,ce1)
        write(*,*)
        write(*,*)'       z           re[e1(z)]        im[e1(z)]'
        write(*,*)' -----------------------------------------------'
        write(*,10)x,y,ce1
10      format(1x,f5.1,2x,f5.1,1x,2d17.8)
        end


        subroutine e1z(z,ce1)
c
c       ====================================================
c       purpose: compute complex exponential integral e1(z)
c       input :  z   --- argument of e1(z)
c       output:  ce1 --- e1(z)
c       ====================================================
c
        implicit complex*16 (c,z)
        implicit double precision (d-h,o-y)
        pi=3.141592653589793d0
        el=0.5772156649015328d0
        x=real(z)
        a0=cdabs(z)
        if (a0.eq.0.0d0) then
           ce1=(1.0d+300,0.0d0)
        else if (a0.le.10.0.or.x.lt.0.0.and.a0.lt.20.0) then
           ce1=(1.0d0,0.0d0)
           cr=(1.0d0,0.0d0)
           do 10 k=1,150
              cr=-cr*k*z/(k+1.0d0)**2
              ce1=ce1+cr
              if (cdabs(cr).le.cdabs(ce1)*1.0d-15) go to 15
10         continue
15         ce1=-el-cdlog(z)+z*ce1
        else
           ct0=(0.0d0,0.0d0)
           do 20 k=120,1,-1
              ct0=k/(1.0d0+k/(z+ct0))
20         continue
           ct=1.0d0/(z+ct0)
           ce1=cdexp(-z)*ct
           if (x.le.0.0.and.dimag(z).eq.0.0) ce1=ce1-pi*(0.0d0,1.0d0)
        endif
        return
        end
