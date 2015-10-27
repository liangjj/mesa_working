        program mcjynb
c
c       ================================================================
c       purpose: this program computes bessel functions jn(z), yn(z)  
c                and their derivatives for a complex argument using
c                subroutine cjynb
c       input :  z --- complex argument of jn(z) and yn(z)
c                n --- order of jn(z) and yn(z)
c                      ( n = 0,1,תתת, n ף 250 )
c       output:  cbj(n) --- jn(z)
c                cdj(n) --- jn'(z)
c                cby(n) --- yn(z)
c                cdy(n) --- yn'(z)
c       eaxmple: z = 4.0 + i 2.0
c
c     n     re[jn(z)]       im[jn(z)]       re[jn'(z)]      im[jn'(z)]
c    -------------------------------------------------------------------
c     0  -.13787022d+01   .39054236d+00   .50735255d+00   .12263041d+01
c     1  -.50735255d+00  -.12263041d+01  -.11546013d+01   .58506793d+00
c     2   .93050039d+00  -.77959350d+00  -.72363400d+00  -.72836666d+00
c     3   .93991546d+00   .23042918d+00   .29742236d+00  -.63587637d+00
c     4   .33565567d+00   .49215925d+00   .47452722d+00  -.29035945d-01
c     5  -.91389835d-02   .28850107d+00   .20054412d+00   .19908868d+00
c
c     n     re[yn(z)]       im[yn(z)]       re[yn'(z)]      im[yn'(z)]
c   --------------------------------------------------------------------
c     0  -.38145893d+00  -.13291649d+01  -.12793101d+01   .51220420d+00
c     1   .12793101d+01  -.51220420d+00  -.58610052d+00  -.10987930d+01
c     2   .79074211d+00   .86842120d+00   .78932897d+00  -.70142425d+00
c     3  -.29934789d+00   .89064431d+00   .70315755d+00   .24423024d+00
c     4  -.61557299d+00   .37996071d+00   .41126221d-01   .34044655d+00
c     5  -.38160033d+00   .20975121d+00  -.33884827d+00  -.20590670d-01
c
c                z = 20.0 + i  10.0 ,      nmax =  5
c
c     n     re[jn(z)]       im[jn(z)]       re[jn'(z)]      im[jn'(z)]
c   --------------------------------------------------------------------
c     0   .15460268d+04  -.10391216d+04  -.10601232d+04  -.15098284d+04
c     1   .10601232d+04   .15098284d+04   .14734253d+04  -.10783122d+04
c     2  -.14008238d+04   .11175029d+04   .11274890d+04   .13643952d+04
c     3  -.11948548d+04  -.12189620d+04  -.11843035d+04   .11920871d+04
c     4   .96778325d+03  -.12666712d+04  -.12483664d+04  -.93887194d+03
c     5   .13018781d+04   .65878188d+03   .64152944d+03  -.12682398d+04
c
c     n     re[yn(z)]       im[yn(z)]       re[yn'(z)]      im[yn'(z)]
c   --------------------------------------------------------------------
c     0   .10391216d+04   .15460268d+04   .15098284d+04  -.10601232d+04
c     1  -.15098284d+04   .10601232d+04   .10783122d+04   .14734253d+04
c     2  -.11175029d+04  -.14008238d+04  -.13643952d+04   .11274890d+04
c     3   .12189620d+04  -.11948548d+04  -.11920871d+04  -.11843035d+04
c     4   .12666712d+04   .96778324d+03   .93887194d+03  -.12483664d+04
c     5  -.65878189d+03   .13018781d+04   .12682398d+04   .64152944d+03
c       ================================================================
c
        implicit double precision (a,b,d-h,o-y)
        implicit complex*16 (c,z)
        common cbj(0:250),cdj(0:250),cby(0:250),cdy(0:250)
        write(*,*)'  please enter n, x,y (z=x+iy) '
        read(*,*)n,x,y
        z=cmplx(x,y)
        write(*,40)x,y,n
        if (n.le.8) then
           ns=1
        else
           write(*,*)'  please enter order step ns'
           read(*,*)ns
        endif
        call cjynb(n,z,nm,cbj,cdj,cby,cdy)
        write(*,*)
        write(*,*)'   n     re[jn(z)]       im[jn(z)]',
     &            '       re[jn''(z)]      im[jn''(z)]'
        write(*,*)' -------------------------------------',
     &            '-------------------------------'
        do 10 k=0,nm,ns
10         write(*,30)k,cbj(k),cdj(k)
        write(*,*)
        write(*,*)'   n     re[yn(z)]       im[yn(z)]',
     &            '       re[yn''(z)]      im[yn''(z)]'
        write(*,*)' -------------------------------------',
     &            '-------------------------------'
        do 20 k=0,nm,ns
20         write(*,30)k,cby(k),cdy(k)
30      format(1x,i4,4d16.8)
40      format(3x,3hz =,f5.1,' + i ',f5.1,' ,',6x,6hnmax =,i3)
        end


        subroutine cjynb(n,z,nm,cbj,cdj,cby,cdy)
c
c       =======================================================
c       purpose: compute bessel functions jn(z), yn(z) and
c                their derivatives for a complex argument
c       input :  z --- complex argument of jn(z) and yn(z)
c                n --- order of jn(z) and yn(z)
c       output:  cbj(n) --- jn(z)
c                cdj(n) --- jn'(z)
c                cby(n) --- yn(z)
c                cdy(n) --- yn'(z)
c                nm --- highest order computed
c       routines called:
c                msta1 and msta2 to calculate the starting
c                point for backward recurrence
c       =======================================================
c
        implicit double precision (a,b,d-h,o-y)
        implicit complex*16 (c,z)
        dimension cbj(0:n),cdj(0:n),cby(0:n),cdy(0:n),
     &            a(4),b(4),a1(4),b1(4)
        el=0.5772156649015329d0
        pi=3.141592653589793d0
        r2p=.63661977236758d0
        y0=dabs(dimag(z))
        a0=cdabs(z)
        nm=n
        if (a0.lt.1.0d-100) then
           do 10 k=0,n
              cbj(k)=(0.0d0,0.0d0)
              cdj(k)=(0.0d0,0.0d0)
              cby(k)=-(1.0d+300,0.0d0)
10            cdy(k)=(1.0d+300,0.0d0)
           cbj(0)=(1.0d0,0.0d0)
           cdj(1)=(0.5d0,0.0d0)
           return
        endif
        if (a0.le.300.d0.or.n.gt.int(0.25*a0)) then
           if (n.eq.0) nm=1
           m=msta1(a0,200)
           if (m.lt.nm) then
              nm=m
           else
              m=msta2(a0,nm,15)
           endif
           cbs=(0.0d0,0.0d0)
           csu=(0.0d0,0.0d0)
           csv=(0.0d0,0.0d0)
           cf2=(0.0d0,0.0d0)
           cf1=(1.0d-100,0.0d0)
           do 15 k=m,0,-1
              cf=2.0d0*(k+1.0d0)/z*cf1-cf2
              if (k.le.nm) cbj(k)=cf
              if (k.eq.2*int(k/2).and.k.ne.0) then
                 if (y0.le.1.0d0) then
                    cbs=cbs+2.0d0*cf
                 else
                    cbs=cbs+(-1)**(k/2)*2.0d0*cf
                 endif
                 csu=csu+(-1)**(k/2)*cf/k
              else if (k.gt.1) then
                 csv=csv+(-1)**(k/2)*k/(k*k-1.0d0)*cf
              endif
              cf2=cf1
15            cf1=cf
           if (y0.le.1.0d0) then
              cs0=cbs+cf
           else
              cs0=(cbs+cf)/cdcos(z)
           endif
           do 20 k=0,nm
20            cbj(k)=cbj(k)/cs0
           ce=cdlog(z/2.0d0)+el
           cby(0)=r2p*(ce*cbj(0)-4.0d0*csu/cs0)
           cby(1)=r2p*(-cbj(0)/z+(ce-1.0d0)*cbj(1)-4.0d0*csv/cs0)
        else
           data a/-.7031250000000000d-01,.1121520996093750d+00,
     &            -.5725014209747314d+00,.6074042001273483d+01/
           data b/ .7324218750000000d-01,-.2271080017089844d+00,
     &             .1727727502584457d+01,-.2438052969955606d+02/
           data a1/.1171875000000000d+00,-.1441955566406250d+00,
     &             .6765925884246826d+00,-.6883914268109947d+01/
           data b1/-.1025390625000000d+00,.2775764465332031d+00,
     &             -.1993531733751297d+01,.2724882731126854d+02/
           ct1=z-0.25d0*pi
           cp0=(1.0d0,0.0d0)
           do 25 k=1,4
25            cp0=cp0+a(k)*z**(-2*k)
           cq0=-0.125d0/z
           do 30 k=1,4
30            cq0=cq0+b(k)*z**(-2*k-1)
           cu=cdsqrt(r2p/z)
           cbj0=cu*(cp0*cdcos(ct1)-cq0*cdsin(ct1))
           cby0=cu*(cp0*cdsin(ct1)+cq0*cdcos(ct1))
           cbj(0)=cbj0
           cby(0)=cby0
           ct2=z-0.75d0*pi
           cp1=(1.0d0,0.0d0)
           do 35 k=1,4
35            cp1=cp1+a1(k)*z**(-2*k)
           cq1=0.375d0/z
           do 40 k=1,4
40            cq1=cq1+b1(k)*z**(-2*k-1)
           cbj1=cu*(cp1*cdcos(ct2)-cq1*cdsin(ct2))
           cby1=cu*(cp1*cdsin(ct2)+cq1*cdcos(ct2))
           cbj(1)=cbj1
           cby(1)=cby1
           do 45 k=2,nm
              cbjk=2.0d0*(k-1.0d0)/z*cbj1-cbj0
              cbj(k)=cbjk
              cbj0=cbj1
45            cbj1=cbjk
        endif
        cdj(0)=-cbj(1)
        do 50 k=1,nm
50         cdj(k)=cbj(k-1)-k/z*cbj(k)
        if (cdabs(cbj(0)).gt.1.0d0) then
           cby(1)=(cbj(1)*cby(0)-2.0d0/(pi*z))/cbj(0)
        endif
        do 55 k=2,nm
           if (cdabs(cbj(k-1)).ge.cdabs(cbj(k-2))) then
              cyy=(cbj(k)*cby(k-1)-2.0d0/(pi*z))/cbj(k-1)
           else
              cyy=(cbj(k)*cby(k-2)-4.0d0*(k-1.0d0)/(pi*z*z))/cbj(k-2)
           endif
           cby(k)=cyy
55      continue
        cdy(0)=-cby(1)
        do 60 k=1,nm
60         cdy(k)=cby(k-1)-k/z*cby(k)
        return
        end


        integer function msta1(x,mp)
c
c       ===================================================
c       purpose: determine the starting point for backward  
c                recurrence such that the magnitude of    
c                jn(x) at that point is about 10^(-mp)
c       input :  x     --- argument of jn(x)
c                mp    --- value of magnitude
c       output:  msta1 --- starting point   
c       ===================================================
c
        implicit double precision (a-h,o-z)
        a0=dabs(x)
        n0=int(1.1*a0)+1
        f0=envj(n0,a0)-mp
        n1=n0+5
        f1=envj(n1,a0)-mp
        do 10 it=1,20             
           nn=n1-(n1-n0)/(1.0d0-f0/f1)                  
           f=envj(nn,a0)-mp
           if(abs(nn-n1).lt.1) go to 20
           n0=n1
           f0=f1
           n1=nn
 10        f1=f
 20     msta1=nn
        return
        end


        integer function msta2(x,n,mp)
c
c       ===================================================
c       purpose: determine the starting point for backward
c                recurrence such that all jn(x) has mp
c                significant digits
c       input :  x  --- argument of jn(x)
c                n  --- order of jn(x)
c                mp --- significant digit
c       output:  msta2 --- starting point
c       ===================================================
c
        implicit double precision (a-h,o-z)
        a0=dabs(x)
        hmp=0.5d0*mp
        ejn=envj(n,a0)
        if (ejn.le.hmp) then
           obj=mp
           n0=int(1.1*a0)
        else
           obj=hmp+ejn
           n0=n
        endif
        f0=envj(n0,a0)-obj
        n1=n0+5
        f1=envj(n1,a0)-obj
        do 10 it=1,20
           nn=n1-(n1-n0)/(1.0d0-f0/f1)
           f=envj(nn,a0)-obj
           if (abs(nn-n1).lt.1) go to 20
           n0=n1
           f0=f1
           n1=nn
10         f1=f
20      msta2=nn+10
        return
        end

        real*8 function envj(n,x)
        double precision x
        envj=0.5d0*dlog10(6.28d0*n)-n*dlog10(1.36d0*x/n)
        return
        end
