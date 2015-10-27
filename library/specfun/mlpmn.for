        program mlpmn
c
c       ==========================================================
c       purpose: this program computes the associated legendre 
c                functions pmn(x) and their derivatives pmn'(x) 
c                using subroutine lpmn
c       input :  x --- argument of pmn(x)
c                m --- order of pmn(x),  m = 0,1,2,...,n
c                n --- degree of pmn(x), n = 0,1,2,...,n
c       output:  pm(m,n) --- pmn(x)
c                pd(m,n) --- pmn'(x)
c       example: x = 0.50
c          pmn(x):
c          m\n        1            2            3            4
c         --------------------------------------------------------
c           0      .500000     -.125000     -.437500     -.289063
c           1     -.866025    -1.299038     -.324760     1.353165
c           2      .000000     2.250000     5.625000     4.218750
c           3      .000000      .000000    -9.742786   -34.099750
c           4      .000000      .000000      .000000    59.062500
c
c          pmn'(x):
c          m\n        1            2            3            4
c         --------------------------------------------------------
c           0     1.000000     1.500000      .375000    -1.562500
c           1      .577350    -1.732051    -6.278684    -5.773503
c           2      .000000    -3.000000     3.750000    33.750000
c           3      .000000      .000000    19.485572      .000000
c           4      .000000      .000000      .000000  -157.500000
c       ==========================================================
c
        implicit double precision (p,x)
        dimension pm(0:100,0:100),pd(0:100,0:100)
        write(*,*)'  please enter m, n and x'
        read(*,*) m,n,x
        write(*,*)
        write(*,*)'  m     n      x          pmn(x)         pmn''(x)'
        write(*,*)' ---------------------------------------------------'
        call lpmn(100,m,n,x,pm,pd)
        do 15 j=0,n
           write(*,10)m,j,x,pm(m,j),pd(m,j)
15      continue
10      format(1x,i3,3x,i3,3x,f5.1,2e17.8)
        end


        subroutine lpmn(mm,m,n,x,pm,pd)
c
c       =====================================================
c       purpose: compute the associated legendre functions 
c                pmn(x) and their derivatives pmn'(x)
c       input :  x  --- argument of pmn(x)
c                m  --- order of pmn(x),  m = 0,1,2,...,n
c                n  --- degree of pmn(x), n = 0,1,2,...,n
c                mm --- physical dimension of pm and pd
c       output:  pm(m,n) --- pmn(x)
c                pd(m,n) --- pmn'(x)
c       =====================================================
c
        implicit double precision (p,x)
        dimension pm(0:mm,0:n),pd(0:mm,0:n)
        do 10 i=0,n
        do 10 j=0,m
           pm(j,i)=0.0d0
10         pd(j,i)=0.0d0
        pm(0,0)=1.0d0
        if (dabs(x).eq.1.0d0) then
           do 15 i=1,n
              pm(0,i)=x**i
15            pd(0,i)=0.5d0*i*(i+1.0d0)*x**(i+1)
           do 20 j=1,n
           do 20 i=1,m
              if (i.eq.1) then
                 pd(i,j)=1.0d+300
              else if (i.eq.2) then
                 pd(i,j)=-0.25d0*(j+2)*(j+1)*j*(j-1)*x**(j+1)
              endif
20         continue
           return
        endif
        ls=1
        if (dabs(x).gt.1.0d0) ls=-1
        xq=dsqrt(ls*(1.0d0-x*x))
        xs=ls*(1.0d0-x*x)
        do 30 i=1,m
30         pm(i,i)=-ls*(2.0d0*i-1.0d0)*xq*pm(i-1,i-1)
        do 35 i=0,m
35         pm(i,i+1)=(2.0d0*i+1.0d0)*x*pm(i,i)
        do 40 i=0,m
        do 40 j=i+2,n
           pm(i,j)=((2.0d0*j-1.0d0)*x*pm(i,j-1)-
     &             (i+j-1.0d0)*pm(i,j-2))/(j-i)
40      continue
        pd(0,0)=0.0d0
        do 45 j=1,n
45         pd(0,j)=ls*j*(pm(0,j-1)-x*pm(0,j))/xs
        do 50 i=1,m
        do 50 j=i,n
           pd(i,j)=ls*i*x*pm(i,j)/xs+(j+i)
     &             *(j-i+1.0d0)/xq*pm(i-1,j)
50      continue
        return
        end
