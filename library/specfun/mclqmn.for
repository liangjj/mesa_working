        program mclqmn
c
c       ============================================================
c       purpose: this program computes the associated legendre 
c                functions qmn(z) and their derivatives qmn'(z) for
c                a complex argument using subroutine clqmn
c       definition: qmn(z)=(-1)**m*(1-z*z)**(m/2)*dm/dzm[qn(z)]
c                   q0(z)=1/2*log[(1+z)/(1-z)]     ( for |z|<1 )
c                   qmn(z)=(z*z-1)**(m/2)*dm/dzm[qn(z)]
c                   q0(z)=1/2*log[(z+1)/(z-1)]     ( for |z|>1 )
c       input :  x --- real part of z
c                y --- imaginary part of z
c                m --- order of qmn(z)  ( m = 0,1,2,תתת )
c                n --- degree of qmn(z) ( n = 0,1,2,תתת )
c       output:  cqm(m,n) --- qmn(z)
c                cqd(m,n) --- qmn'(z)
c       examples:
c                n = 5, x = 0.5, y = 0.2
c
c       m     re[qmn(z)]    im[qmn(z)]    re[qmn'(z)]   im[qmn'(z)]
c      -------------------------------------------------------------
c       0    .987156d+00   .354345d+00    .324023d+01  -.447297d+01
c       1   -.240328d+01   .436861d+01    .281158d+02   .171437d+02
c       2   -.245853d+02  -.138072d+02   -.106283d+03   .913792d+02
c       3    .102723d+03  -.651233d+02   -.362578d+03  -.429802d+03
c       4    .155510d+03   .357712d+03    .196975d+04  -.287414d+02
c       5   -.167357d+04  -.680954d+03   -.193093d+04  -.925757d+03
c
c                n = 5, x = 2.5, y = 1.0
c
c       m     re[qmn(z)]    im[qmn(z)]    re[qmn'(z)]   im[qmn'(z)]
c      -------------------------------------------------------------
c       0   -.274023d-04  -.227141d-04    .809834d-04   .210884d-04
c       1    .165620d-03   .136108d-03   -.489095d-03  -.124400d-03
c       2   -.118481d-02  -.948832d-03    .349090d-02   .825057d-03
c       3    .982179d-02   .753264d-02   -.288271d-01  -.596384d-02
c       4   -.927915d-01  -.669521d-01    .270840d+00   .451376d-01
c       5    .985601d+00   .656737d+00   -.285567d+01  -.332533d+00
c       ============================================================
c
        implicit double precision (x,y)
        implicit complex*16 (c,z)
        dimension cqm(0:40,0:40),cqd(0:40,0:40)
        write(*,*)'  please enter m, n, x and y '
        read(*,*) m,n,x,y
        write(*,30)m,n,x,y
        call clqmn(40,m,n,x,y,cqm,cqd)
        write(*,*)'   m   n   re[qmn(z)]    im[qmn(z)]    ',
     &            're[qmn''(z)]   im[qmn''(z)]'
        write(*,*)' -----------------------------------',
     &            '------------------------------'
        do 10 j=0,n
10         write(*,20)m,j,cqm(m,j),cqd(m,j)
20      format(1x,2i4,2d14.6,1x,2d14.6)
30      format(1x,'m =',i2,', ','n =',i2,', ','x =',f4.1,
     &         ', ','y =',f4.1)
        end


        subroutine clqmn(mm,m,n,x,y,cqm,cqd)
c
c       =======================================================
c       purpose: compute the associated legendre functions of
c                the second kind, qmn(z) and qmn'(z), for a
c                complex argument
c       input :  x  --- real part of z
c                y  --- imaginary part of z
c                m  --- order of qmn(z)  ( m = 0,1,2,תתת )
c                n  --- degree of qmn(z) ( n = 0,1,2,תתת )
c                mm --- physical dimension of cqm and cqd
c       output:  cqm(m,n) --- qmn(z)
c                cqd(m,n) --- qmn'(z)
c       =======================================================
c
        implicit double precision (x,y)
        implicit complex*16 (c,z)
        dimension cqm(0:mm,0:n),cqd(0:mm,0:n)
        z=cmplx(x,y)
        if (dabs(x).eq.1.0d0.and.y.eq.0.0d0) then
           do 10 i=0,m
           do 10 j=0,n
              cqm(i,j)=(1.0d+300,0.0d0)
              cqd(i,j)=(1.0d+300,0.0d0)
10         continue
           return
        endif
        xc=cdabs(z)
        if (dimag(z).eq.0.0d0.or.xc.lt.1.0d0) ls=1
        if (xc.gt.1.0d0) ls=-1
        zq=cdsqrt(ls*(1.0d0-z*z))
        zs=ls*(1.0d0-z*z)
        cq0=0.5d0*cdlog(ls*(1.0d0+z)/(1.0d0-z))
        if (xc.lt.1.0001d0) then
           cqm(0,0)=cq0
           cqm(0,1)=z*cq0-1.0d0
           cqm(1,0)=-1.0d0/zq
           cqm(1,1)=-zq*(cq0+z/(1.0d0-z*z))
           do 15 i=0,1
           do 15 j=2,n
              cqm(i,j)=((2.0d0*j-1.0d0)*z*cqm(i,j-1)
     &                -(j+i-1.0d0)*cqm(i,j-2))/(j-i)
15         continue
           do 20 j=0,n
           do 20 i=2,m
              cqm(i,j)=-2.0d0*(i-1.0d0)*z/zq*cqm(i-1,j)-ls*
     &                 (j+i-1.0d0)*(j-i+2.0d0)*cqm(i-2,j)
20         continue
        else
           if (xc.gt.1.1) then
              km=40+m+n
           else
              km=(40+m+n)*int(-1.0-1.8*log(xc-1.0))
           endif
           cqf2=(0.0d0,0.0d0)
           cqf1=(1.0d0,0.0d0)
           do 25 k=km,0,-1
              cqf0=((2*k+3.0d0)*z*cqf1-(k+2.0d0)*cqf2)/(k+1.0d0)
              if (k.le.n) cqm(0,k)=cqf0
              cqf2=cqf1
25            cqf1=cqf0
           do 30 k=0,n
30            cqm(0,k)=cq0*cqm(0,k)/cqf0
           cqf2=0.0d0
           cqf1=1.0d0
           do 35 k=km,0,-1
              cqf0=((2*k+3.0d0)*z*cqf1-(k+1.0d0)*cqf2)/(k+2.0d0)
              if (k.le.n) cqm(1,k)=cqf0
              cqf2=cqf1
35            cqf1=cqf0
           cq10=-1.0d0/zq
           do 40 k=0,n
40            cqm(1,k)=cq10*cqm(1,k)/cqf0
           do 45 j=0,n
              cq0=cqm(0,j)
              cq1=cqm(1,j)
              do 45 i=0,m-2
                 cqf=-2.0d0*(i+1)*z/zq*cq1+(j-i)*(j+i+1.0d0)*cq0
                 cqm(i+2,j)=cqf
                 cq0=cq1
                 cq1=cqf
45         continue
        endif
        cqd(0,0)=ls/zs
        do 50 j=1,n
50         cqd(0,j)=ls*j*(cqm(0,j-1)-z*cqm(0,j))/zs
        do 55 j=0,n
        do 55 i=1,m
           cqd(i,j)=ls*i*z/zs*cqm(i,j)+(i+j)*(j-i+1.0d0)
     &              /zq*cqm(i-1,j)
55      continue
        return
        end


