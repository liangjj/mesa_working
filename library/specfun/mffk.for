        program mffk
c
c       ==============================================================
c       purpose: this program computes the modified fresnel integrals 
c                fñ(x) and kñ(x) using subroutine ffk
c       input :  x   --- argument of fñ(x) and kñ(x)
c                ks  --- sign code
c                        ks=0 for calculating f+(x) and k+(x)
c                        ks=1 for calculating f_(x) and k_(x)
c       output:  fr  --- re[fñ(x)]
c                fi  --- im[fñ(x)]
c                fm  --- |fñ(x)|
c                fa  --- arg[fñ(x)]  (degs.)
c                gr  --- re[kñ(x)]
c                gi  --- im[kñ(x)]
c                gm  --- |kñ(x)|
c                ga  --- arg[kñ(x)]  (degs.)
c       example:
c
c         x     re[fñ(x)]   ñim[fñ(x)]    mod[fñ(x)]  ñarg[fñ(x)]
c       ----------------------------------------------------------
c        0.0    .62665707    .62665707    .88622693    45.000000
c        2.0    .16519561   -.17811942    .24293233   -47.155835
c        4.0    .03219674   -.12047678    .12470479   -75.037684
c        6.0    .08245304   -.01180212    .08329342    -8.145843
c        8.0   -.05729996    .02493542    .06249048   156.482601
c       10.0    .02553188    .04298617    .04999688    59.291561
c
c         x     re[kñ(x)]   ñim[kñ(x)]    mod[kñ(x)]  ñarg[kñ(x)]
c       ----------------------------------------------------------
c        0.0    .50000000    .00000000    .50000000     0.000000
c        2.0    .10702394    .08562295    .13705989    38.661047
c        4.0    .05126306    .04818949    .07035714    43.229843
c        6.0    .03368650    .03276566    .04699328    44.206095
c        8.0    .02512396    .02473472    .03525648    44.552712
c       10.0    .02004532    .01984592    .02820772    44.713609
c       ===============================================================
c
        implicit double precision (f,g,x)
        write(*,*)'please enter x'
        read(*,*) x
        write(*,*)
        write(*,*)'   x      re[fñ(x)]    ñim[fñ(x)]     ',
     &             'mod[fñ(x)]   ñarg[fñ(x)]'
        write(*,*)' ---------------------------------------',
     &            '-----------------------'
        call ffk(0,x,fr,fi,fm,fa,gr,gi,gm,ga)
        write(*,10)x,fr,fi,fm,fa
        write(*,*)
        write(*,*)'   x      re[kñ(x)]    ñim[kñ(x)]     ',
     &             'mod[kñ(x)]   ñarg[kñ(x)]'
        write(*,*)' ---------------------------------------',
     &            '-----------------------'
        write(*,10)x,gr,gi,gm,ga
10      format(1x,f5.1,3f14.8,f14.6)
        end


        subroutine ffk(ks,x,fr,fi,fm,fa,gr,gi,gm,ga)
c
c       =======================================================
c       purpose: compute modified fresnel integrals fñ(x) 
c                and kñ(x)
c       input :  x   --- argument of fñ(x) and kñ(x)
c                ks  --- sign code
c                        ks=0 for calculating f+(x) and k+(x)
c                        ks=1 for calculating f_(x) and k_(x)
c       output:  fr  --- re[fñ(x)]
c                fi  --- im[fñ(x)]
c                fm  --- |fñ(x)|
c                fa  --- arg[fñ(x)]  (degs.)
c                gr  --- re[kñ(x)]
c                gi  --- im[kñ(x)]
c                gm  --- |kñ(x)|
c                ga  --- arg[kñ(x)]  (degs.)
c       ======================================================
c
        implicit double precision (a-h,o-z)
        srd= 57.29577951308233d0
        eps=1.0d-15
        pi=3.141592653589793d0
        pp2=1.2533141373155d0
        p2p=.7978845608028654d0
        xa=dabs(x)
        x2=x*x
        x4=x2*x2
        if (x.eq.0.0d0) then
           fr=.5d0*dsqrt(0.5d0*pi)
           fi=(-1)**ks*fr
           fm=dsqrt(0.25d0*pi)
           fa=(-1)**ks*45.0d0
           gr=.5d0
           gi=0.0d0
           gm=.5d0
           ga=0.0d0
        else
           if (xa.le.2.5d0) then
              xr=p2p*xa
              c1=xr
              do 10 k=1,50
                 xr=-.5d0*xr*(4.0d0*k-3.0d0)/k/(2.0d0*k-1.0d0)
     &              /(4.0d0*k+1.0d0)*x4
                 c1=c1+xr
                 if (dabs(xr/c1).lt.eps) go to 15
10            continue
15            s1=p2p*xa*xa*xa/3.0d0
              xr=s1
              do 20 k=1,50
                 xr=-.5d0*xr*(4.0d0*k-1.0d0)/k/(2.0d0*k+1.0d0)
     &              /(4.0d0*k+3.0d0)*x4
                 s1=s1+xr
                 if (dabs(xr/s1).lt.eps) go to 40
20            continue
           else if (xa.lt.5.5d0) then
              m=int(42+1.75*x2)
              xsu=0.0d0
              xc=0.0d0
              xs=0.0d0
              xf1=0.0d0
              xf0=1d-100
              do 25 k=m,0,-1
                 xf=(2.0d0*k+3.0d0)*xf0/x2-xf1
                 if (k.eq.2*int(k/2))  then
                    xc=xc+xf
                 else
                    xs=xs+xf
                 endif
                 xsu=xsu+(2.0d0*k+1.0d0)*xf*xf
                 xf1=xf0
25               xf0=xf
              xq=dsqrt(xsu)
              xw=p2p*xa/xq
              c1=xc*xw
              s1=xs*xw
           else
              xr=1.0d0
              xf=1.0d0
              do 30 k=1,12
                 xr=-.25d0*xr*(4.0d0*k-1.0d0)*(4.0d0*k-3.0d0)/x4
30               xf=xf+xr
              xr=1.0d0/(2.0d0*xa*xa)
              xg=xr
              do 35 k=1,12
                 xr=-.25d0*xr*(4.0d0*k+1.0d0)*(4.0d0*k-1.0d0)/x4
35               xg=xg+xr
              c1=.5d0+(xf*dsin(x2)-xg*dcos(x2))/dsqrt(2.0d0*pi)/xa
              s1=.5d0-(xf*dcos(x2)+xg*dsin(x2))/dsqrt(2.0d0*pi)/xa
           endif
40         fr=pp2*(.5d0-c1)
           fi0=pp2*(.5d0-s1)
           fi=(-1)**ks*fi0
           fm=dsqrt(fr*fr+fi*fi)
           if (fr.ge.0.0) then
              fa=srd*datan(fi/fr)
           else if (fi.gt.0.0) then
              fa=srd*(datan(fi/fr)+pi)
           else if (fi.lt.0.0) then
              fa=srd*(datan(fi/fr)-pi)
           endif
           xp=x*x+pi/4.0d0
           cs=dcos(xp)
           ss=dsin(xp)
           xq2=1.0d0/dsqrt(pi)
           gr=xq2*(fr*cs+fi0*ss)
           gi=(-1)**ks*xq2*(fi0*cs-fr*ss)
           gm=dsqrt(gr*gr+gi*gi)
           if (gr.ge.0.0) then
              ga=srd*datan(gi/gr)
           else if (gi.gt.0.0) then
              ga=srd*(datan(gi/gr)+pi)
           else if (gi.lt.0.0) then
              ga=srd*(datan(gi/gr)-pi)
           endif
           if (x.lt.0.0d0) then
              fr=pp2-fr
              fi=(-1)**ks*pp2-fi
              fm=dsqrt(fr*fr+fi*fi)
              fa=srd*datan(fi/fr)
              gr=dcos(x*x)-gr
              gi=-(-1)**ks*dsin(x*x)-gi
              gm=dsqrt(gr*gr+gi*gi)
              ga=srd*datan(gi/gr)
           endif
        endif
        return
        end
