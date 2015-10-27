        program mclqn
c
c       ==========================================================
c       purpose: this program computes the legendre polynomials 
c                qn(z) and qn'(z) for a complex argument using
c                subroutine clqn
c       input :  x --- real part of z
c                y --- imaginary part of z
c                n --- degree of qn(z), n = 0,1,...
c       output:  cqn(n) --- qn(z)
c                cqd(n) --- qn'(z)
c       examples:
c
c       z = 0.5 + 0.5 i
c       n    re[qn(z)]     im[qn(z)]     re[qn'(z)]    im[qn'(z)]
c      -----------------------------------------------------------
c       0   .402359d+00   .553574d+00   .800000d+00   .400000d+00
c       1  -.107561d+01   .477967d+00   .602359d+00   .115357d+01
c       2  -.136636d+01  -.725018d+00  -.242682d+01   .183390d+01
c       3   .182619d+00  -.206146d+01  -.622944d+01  -.247151d+01
c       4   .298834d+01  -.110022d+01  -.114849d+01  -.125963d+02
c       5   .353361d+01   .334847d+01   .206656d+02  -.123735d+02
c
c       z = 3.0 + 2.0 i
c       n    re[qn(z)]     im[qn(z)]     re[qn'(z)]    im[qn'(z)]
c      -----------------------------------------------------------
c       0   .229073d+00  -.160875d+00  -.250000d-01   .750000d-01
c       1   .896860d-02  -.244805d-01   .407268d-02   .141247d-01
c       2  -.736230d-03  -.281865d-02   .190581d-02   .155860d-02
c       3  -.264727d-03  -.227023d-03   .391535d-03   .314880d-04
c       4  -.430648d-04  -.443187d-05   .527190d-04  -.305592d-04
c       5  -.481362d-05   .265297d-05   .395108d-05  -.839883d-05
c       ==========================================================
c
        implicit double precision (x,y)
        implicit complex*16 (c,z)
        dimension cqn(0:100),cqd(0:100)
        write(*,*)'  please enter nmax, x and y (z=x+iy)'
        read(*,*)n,x,y
        write(*,30)x,y
        write(*,*)
        call clqn(n,x,y,cqn,cqd)
        write(*,*)'  n    re[qn(z)]     im[qn(z)]     re[qn''(z)]',
     &            '    im[qn''(z)]'
        write(*,*)' ---------------------------------------------',
     &            '--------------'
        do 10 k=0,n
10         write(*,20)k,cqn(k),cqd(k)
20      format(1x,i3,4d14.6)
30      format(3x,'x =',f5.1,',  ','y =',f5.1)
        end


        subroutine clqn(n,x,y,cqn,cqd)
c
c       ==================================================
c       purpose: compute the legendre functions qn(z) and
c                their derivatives qn'(z) for a complex
c                argument
c       input :  x --- real part of z
c                y --- imaginary part of z
c                n --- degree of qn(z), n = 0,1,2,...
c       output:  cqn(n) --- qn(z)
c                cqd(n) --- qn'(z)
c       ==================================================
c
        implicit double precision (x,y)
        implicit complex*16 (c,z)
        dimension cqn(0:n),cqd(0:n)
        z=cmplx(x,y)
        if (z.eq.1.0d0) then
           do 10 k=0,n
              cqn(k)=(1.0d+300,0.0d0)
10            cqd(k)=(1.0d+300,0.0d0)
           return
        endif
        ls=1
        if (cdabs(z).gt.1.0d0) ls=-1
        cq0=0.5d0*cdlog(ls*(1.0d0+z)/(1.0d0-z))
        cq1=z*cq0-1.0d0
        cqn(0)=cq0
        cqn(1)=cq1
        if (cdabs(z).lt.1.0001d0) then
           cqf0=cq0
           cqf1=cq1
           do 15 k=2,n
              cqf2=((2.0d0*k-1.0d0)*z*cqf1-(k-1.0d0)*cqf0)/k
              cqn(k)=cqf2
              cqf0=cqf1
15            cqf1=cqf2
        else
           if (cdabs(z).gt.1.1d0) then
              km=40+n
           else
              km=(40+n)*int(-1.0-1.8*log(cdabs(z-1.0)))
           endif
           cqf2=0.0d0
           cqf1=1.0d0
           do 20 k=km,0,-1
              cqf0=((2*k+3.0d0)*z*cqf1-(k+2.0d0)*cqf2)/(k+1.0d0)
              if (k.le.n) cqn(k)=cqf0
              cqf2=cqf1
20            cqf1=cqf0
           do 25 k=0,n
25            cqn(k)=cqn(k)*cq0/cqf0
        endif
        cqd(0)=(cqn(1)-z*cqn(0))/(z*z-1.0d0)
        do 30 k=1,n
30         cqd(k)=(k*z*cqn(k)-k*cqn(k-1))/(z*z-1.0d0)
        return
        end
