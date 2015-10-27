        program mjdzo
c
c       =============================================================
c       purpose: this program computes the zeros of bessel functions 
c                jn(x) and jn'(x), and arranges them in the order 
c                of their values 
c       input :  nt    --- number of total zeros ( nt ף 1200 )
c       output:  zo(l) --- value of the l-th zero of jn(x) and 
c                          jn'(x)
c                n(l)  --- n, order of jn(x) or jn'(x) associated
c                          with the l-th zero
c                m(l)  --- m, serial number of the zeros of jn(x)
c                          or jn'(x) associated with the l-th zero
c                          ( l is the serial number of all the
c                            zeros of jn(x) and jn'(x) )
c                p(l)  --- tm or te, a code for designating the
c                          zeros of jn(x) or jn'(x)
c                          in the waveguide applications, the zeros
c                          of jn(x) correspond to tm modes and those
c                          of jn'(x) correspond to te modes.
c       =============================================================
c
        implicit double precision (a-h,o-z)
        character p(1400)*4
        dimension n(1400),m(1400),zo(1400)
        write(*,*)'nt=?'
        read(*,*)nt
        write(*,60)nt
        write(*,70)
        call jdzo(nt,n,m,p,zo)
        write(*,*)
        ks=nt/101+1
        do 55 k0=1,ks
           write(*,*)' table           zeros of bessel',
     &               ' functions jn(x) and jn''(x)'
           write(*,*)
           write(*,*)' ----------------------------------',
     &               '----------------------------------'
           do 50 k=1,50
              j1=100*(k0-1)+k+1
              j2=j1+50
              if (j1.le.nt+1.and.j2.le.nt+1) then
                  write(*,65)j1-1,p(j1),n(j1),m(j1),zo(j1),
     &                 j2-1,p(j2),n(j2),m(j2),zo(j2)
              else if (j1.le.nt+1.and.j2.gt.nt+1) then
                  write(*,65)j1-1,p(j1),n(j1),m(j1),zo(j1)
              endif
50         continue
           write(*,*)' ----------------------------------',
     &               '----------------------------------'
           write(*,75)
55      continue
60      format(1x,'total number of the zeros:',i5)
65      format(1x,i4,3x,a2,i4,2h -,i2,f14.8,3x,1h|,2x,i4,
     &         3x,a2,i4,2h -,i2,f14.8)
70      format(15x,'***  please wait !  the program is running  ***')
75      format(/)
        end


        subroutine jdzo(nt,n,m,p,zo)
c
c       ===========================================================
c       purpose: compute the zeros of bessel functions jn(x) and
c                jn'(x), and arrange them in the order of their
c                magnitudes
c       input :  nt    --- number of total zeros ( nt ף 1200 )
c       output:  zo(l) --- value of the l-th zero of jn(x)
c                          and jn'(x)
c                n(l)  --- n, order of jn(x) or jn'(x) associated
c                          with the l-th zero
c                m(l)  --- m, serial number of the zeros of jn(x)
c                          or jn'(x) associated with the l-th zero
c                          ( l is the serial number of all the
c                            zeros of jn(x) and jn'(x) )
c                p(l)  --- tm or te, a code for designating the
c                          zeros of jn(x)  or jn'(x).
c                          in the waveguide applications, the zeros
c                          of jn(x) correspond to tm modes and 
c                          those of jn'(x) correspond to te modes
c       routine called:    bjndd for computing jn(x), jn'(x) and  
c                          jn''(x)
c       =============================================================
c
        implicit double precision (a-h,o-z)
        character p(1400)*4,p1(70)*4
        dimension n(1400),m(1400),zo(1400),n1(70),m1(70),
     &            zoc(70),bj(101),dj(101),fj(101)
        if (nt.lt.600) then
           xm=-1.0+2.248485*nt**0.5-.0159382*nt+3.208775e-4
     &        *nt**1.5
           nm=int(14.5+.05875*nt)
           mm=int(.02*nt)+6
        else
           xm=5.0+1.445389*nt**.5+.01889876*nt-2.147763e-4
     &        *nt**1.5
           nm=int(27.8+.0327*nt)
           mm=int(.01088*nt)+10
        endif
        l0=0
        do 45 i=1,nm
           x1=.407658+.4795504*(i-1)**.5+.983618*(i-1)
           x2=1.99535+.8333883*(i-1)**.5+.984584*(i-1)
           l1=0
           do 30 j=1,mm
              if (i.eq.1.and.j.eq.1) go to 15
              x=x1
10            call bjndd(i,x,bj,dj,fj)
              x0=x
              x=x-dj(i)/fj(i)
              if (x1.gt.xm) go to 20
              if (dabs(x-x0).gt.1.0d-10) go to 10
15            l1=l1+1
              n1(l1)=i-1
              m1(l1)=j
              if (i.eq.1) m1(l1)=j-1
              p1(l1)='te'
              zoc(l1)=x
              if (i.le.15) then
                 x1=x+3.057+.0122*(i-1)+(1.555+.41575*(i-1))/(j+1)**2
              else
                 x1=x+2.918+.01924*(i-1)+(6.26+.13205*(i-1))/(j+1)**2
              endif
20            x=x2
25            call bjndd(i,x,bj,dj,fj)
              x0=x
              x=x-bj(i)/dj(i)
              if (x.gt.xm) go to 30
              if (dabs(x-x0).gt.1.0d-10) go to 25
              l1=l1+1
              n1(l1)=i-1
              m1(l1)=j
              p1(l1)='tm'
              zoc(l1)=x
              if (i.le.15) then
                 x2=x+3.11+.0138*(i-1)+(.04832+.2804*(i-1))/(j+1)**2
              else
                 x2=x+3.001+.0105*(i-1)+(11.52+.48525*(i-1))/(j+3)**2
              endif
30         continue
           l=l0+l1
           l2=l
35         if (l0.eq.0) then
              do 40 k=1,l
                 zo(k)=zoc(k)
                 n(k)=n1(k)
                 m(k)=m1(k)
40               p(k)=p1(k)
              l1=0
           else if (l0.ne.0) then
              if (zo(l0).ge.zoc(l1)) then
                 zo(l0+l1)=zo(l0)
                 n(l0+l1)=n(l0)
                 m(l0+l1)=m(l0)
                 p(l0+l1)=p(l0)
                 l0=l0-1
              else
                 zo(l0+l1)=zoc(l1)
                 n(l0+l1)=n1(l1)
                 m(l0+l1)=m1(l1)
                 p(l0+l1)=p1(l1)
                 l1=l1-1
              endif
           endif
           if (l1.ne.0) go to 35
45         l0=l2
        return
        end


        subroutine bjndd(n,x,bj,dj,fj)
c
c       =====================================================
c       purpose: compute bessel functions jn(x) and their
c                first and second derivatives ( n= 0,1,תתת )
c       input:   x ---  argument of jn(x)  ( x ע 0 )
c                n ---  order of jn(x)
c       output:  bj(n+1) ---  jn(x)
c                dj(n+1) ---  jn'(x)
c                fj(n+1) ---  jn"(x)
c       =====================================================
c
        implicit double precision (a-h,o-z)
        dimension bj(101),dj(101),fj(101)
        do 10 nt=1,900
           mt=int(0.5*log10(6.28*nt)-nt*log10(1.36*dabs(x)/nt))
           if (mt.gt.20) go to 15
10      continue
15      m=nt
        bs=0.0d0
        f0=0.0d0
        f1=1.0d-35
        do 20 k=m,0,-1
           f=2.0d0*(k+1.0d0)*f1/x-f0
           if (k.le.n) bj(k+1)=f
           if (k.eq.2*int(k/2)) bs=bs+2.0d0*f
           f0=f1
20         f1=f
        do 25 k=0,n
25         bj(k+1)=bj(k+1)/(bs-f)
        dj(1)=-bj(2)
        fj(1)=-1.0d0*bj(1)-dj(1)/x
        do 30 k=1,n
           dj(k+1)=bj(k)-k*bj(k+1)/x
30         fj(k+1)=(k*k/(x*x)-1.0d0)*bj(k+1)-dj(k+1)/x
        return
        end
