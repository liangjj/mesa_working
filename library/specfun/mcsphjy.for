        program mcsphjy
c
c       ================================================================
c       purpose: this program computes the spherical bessel functions 
c                jn(z), yn(z), and their derivatives for a complex
c                argument using subroutine csphjy
c       input :  z --- complex argument
c                n --- order of jn(z) & yn(z) ( 0 ó n ó 250 )
c       output:  csj(n) --- jn(z)
c                cdj(n) --- jn'(z)
c                csy(n) --- yn(z)
c                cdy(n) --- yn'(z)
c       example: z = 4.0+i 2.0
c
c     n     re[jn(z)]       im[jn(z)]       re[jn'(z)]      im[jn'(z)]
c   --------------------------------------------------------------------
c     0  -.80651523d+00  -.18941093d+00  -.37101203d-01   .75210758d+00
c     1   .37101203d-01  -.75210758d+00  -.67093420d+00   .11885235d+00
c     2   .60314368d+00  -.27298399d+00  -.24288981d+00  -.40737409d+00
c     3   .42955048d+00   .17755176d+00   .18848259d+00  -.24320520d+00
c     4   .12251323d+00   .22087111d+00   .19660170d+00   .17937264d-01
c     5  -.10242676d-01   .10975433d+00   .68951842d-01   .83020305d-01
c
c     n     re[yn(z)]       im[yn(z)]       re[yn'(z)]      im[yn'(z)]
c   --------------------------------------------------------------------
c     0   .21734534d+00  -.79487692d+00  -.77049661d+00  -.87010064d-02
c     1   .77049661d+00   .87010064d-02  -.92593503d-01  -.64425800d+00
c     2   .24756293d+00   .56894854d+00   .45127429d+00  -.25839924d+00
c     3  -.23845941d+00   .43646607d+00   .26374403d+00   .12439192d+00
c     4  -.27587985d+00   .20902555d+00  -.67092335d-01   .89500599d-01
c     5  -.70001327d-01   .18807178d+00  -.30472133d+00  -.58661384d-01
c       ================================================================
c
        implicit complex*16 (c,z)
        double precision x,y
        dimension csj(0:250),cdj(0:250),csy(0:250),cdy(0:250)
        write(*,*)'please enter n,x,y (z=x+iy) '
        read(*,*)n,x,y
        write(*,30)n,x,y
        z=cmplx(x,y)
        if (n.le.8) then
           ns=1
        else
           write(*,*)'please enter order step ns '
           read(*,*)ns
        endif
        call csphjy(n,z,nm,csj,cdj,csy,cdy)
        write(*,*)
        write(*,*)'  n      re[jn(z)]        im[jn(z)]',
     &  '        re[jn''(z)]       im[jn''(z)]'
        write(*,*)'--------------------------------------------',
     &  '----------------------------'
        do 10 k=0,nm,ns
10         write(*,20)k,csj(k),cdj(k)
        write(*,*)
        write(*,*)'  n      re[yn(z)]        im[yn(z)]',
     &  '        re[yn''(z)]       im[yn''(z)]'
        write(*,*)'--------------------------------------------',
     &  '----------------------------'
        do 15 k=0,nm,ns
15         write(*,20)k,csy(k),cdy(k)
20      format(1x,i3,4d17.8)
30      format(3x,6hnmaz =,i3,',     ','z = ',f8.1,'+ i',f8.1)
        end


        subroutine csphjy(n,z,nm,csj,cdj,csy,cdy)
c
c       ==========================================================
c       purpose: compute spherical bessel functions jn(z) & yn(z)
c                and their derivatives for a complex argument
c       input :  z --- complex argument
c                n --- order of jn(z) & yn(z) ( n = 0,1,2,... )
c       output:  csj(n) --- jn(z)
c                cdj(n) --- jn'(z)
c                csy(n) --- yn(z)
c                cdy(n) --- yn'(z)
c                nm --- highest order computed
c       routines called:
c                msta1 and msta2 for computing the starting
c                point for backward recurrence
c       ==========================================================
c
        implicit complex*16 (c,z)
        double precision a0
        dimension csj(0:n),cdj(0:n),csy(0:n),cdy(0:n)
        a0=cdabs(z)
        nm=n
        if (a0.lt.1.0d-60) then
           do 10 k=0,n
              csj(k)=0.0d0
              cdj(k)=0.0d0
              csy(k)=-1.0d+300
10            cdy(k)=1.0d+300
           csj(0)=(1.0d0,0.0d0)
           cdj(1)=(.333333333333333d0,0.0d0)
           return
        endif
        csj(0)=cdsin(z)/z
        csj(1)=(csj(0)-cdcos(z))/z
        if (n.ge.2) then
           csa=csj(0)
           csb=csj(1)
           m=msta1(a0,200)
           if (m.lt.n) then
              nm=m
           else
              m=msta2(a0,n,15)
           endif
           cf0=0.0d0
           cf1=1.0d0-100
           do 15 k=m,0,-1
              cf=(2.0d0*k+3.0d0)*cf1/z-cf0
              if (k.le.nm) csj(k)=cf
              cf0=cf1
15            cf1=cf
           if (cdabs(csa).gt.cdabs(csb)) cs=csa/cf
           if (cdabs(csa).le.cdabs(csb)) cs=csb/cf0
           do 20 k=0,nm
20            csj(k)=cs*csj(k)
        endif
        cdj(0)=(cdcos(z)-cdsin(z)/z)/z
        do 25 k=1,nm
25         cdj(k)=csj(k-1)-(k+1.0d0)*csj(k)/z
        csy(0)=-cdcos(z)/z
        csy(1)=(csy(0)-cdsin(z))/z
        cdy(0)=(cdsin(z)+cdcos(z)/z)/z
        cdy(1)=(2.0d0*cdy(0)-cdcos(z))/z
        do 30 k=2,nm
           if (cdabs(csj(k-1)).gt.cdabs(csj(k-2))) then
              csy(k)=(csj(k)*csy(k-1)-1.0d0/(z*z))/csj(k-1)
           else
              csy(k)=(csj(k)*csy(k-2)-(2.0d0*k-1.0d0)/z**3)/csj(k-2)
           endif
30      continue
        do 35 k=2,nm
35         cdy(k)=csy(k-1)-(k+1.0d0)*csy(k)/z
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
