        program mfcszo
c
c       ===========================================================
c       purpose : this program computes the complex zeros of the
c                 fresnel integral c(z) or s(z) using subroutine
c                 fcszo
c       input :   kf  --- function code
c                         kf=1 for c(z) or kf=2 for s(z)
c                 nt  --- total number of zeros
c       output:   zo(l) --- l-th zero of c(z) or s(z)
c       example:  nt=10
c
c       n     complex zeros of c(z)        complex zeros of s(z)
c     ------------------------------------------------------------
c       1    1.7436675 + i .30573506      2.0092570 + i .28854790
c       2    2.6514596 + i .25290396      2.8334772 + i .24428524
c       3    3.3203593 + i .22395346      3.4675331 + i .21849268
c       4    3.8757345 + i .20474747      4.0025782 + i .20085103
c       5    4.3610635 + i .19066973      4.4741893 + i .18768859
c       6    4.7976077 + i .17970801      4.9006784 + i .17732036
c       7    5.1976532 + i .17081930      5.2929467 + i .16884418
c       8    5.5690602 + i .16339854      5.6581068 + i .16172492
c       9    5.9172173 + i .15706585      6.0011034 + i .15562108
c      10    6.2460098 + i .15156826      6.3255396 + i .15030246
c       ===========================================================
c
        implicit double precision (e,p,w)
        implicit complex *16 (c,z)
        dimension zo(100)
        write(*,*)'please enter kf and nt '
        read(*,*)kf,nt
        write(*,20)kf,nt
        write(*,*)' *****     please wait !     *****'
        call fcszo(kf,nt,zo)
        write(*,*)
        if (kf.eq.1) write(*,*)'  n        complex zeros of c(z)'
        if (kf.eq.2) write(*,*)'  n        complex zeros of s(z)'
        write(*,*)'-----------------------------------'
        do 10 i=1,nt
10         write(*,30) i,zo(i)
20      format(2x,'kf=',i2,',     ','nt=',i3)
30      format(1x,i3,f13.8,2x,2h+i,f13.8)
        end


        subroutine fcszo(kf,nt,zo)
c
c       ===============================================================
c       purpose: compute the complex zeros of fresnel integral c(z)
c                or s(z) using modified newton's iteration method
c       input :  kf  --- function code
c                        kf=1 for c(z) or kf=2 for s(z)
c                nt  --- total number of zeros
c       output:  zo(l) --- l-th zero of c(z) or s(z)
c       routines called: 
c            (1) cfc for computing fresnel integral c(z)
c            (2) cfs for computing fresnel integral s(z)
c       ==============================================================
c
        implicit double precision (e,p,w)
        implicit complex *16 (c,z)
        dimension zo(nt)
        pi=3.141592653589793d0
        do 35 nr=1,nt
           if (kf.eq.1) psq=dsqrt(4.0d0*nr-1.0d0)
           if (kf.eq.2) psq=2.0d0*nr**(0.5)
           px=psq-dlog(pi*psq)/(pi*pi*psq**3.0)
           py=dlog(pi*psq)/(pi*psq)
           z=cmplx(px,py)
           if (kf.eq.2) then
              if (nr.eq.2) z=(2.8334,0.2443)
              if (nr.eq.3) z=(3.4674,0.2185)
              if (nr.eq.4) z=(4.0025,0.2008)
           endif
           it=0
15         it=it+1
           if (kf.eq.1) call cfc(z,zf,zd)
           if (kf.eq.2) call cfs(z,zf,zd)
           zp=(1.0d0,0.0d0)
           do 20 i=1,nr-1
20            zp=zp*(z-zo(i))
           zfd=zf/zp
           zq=(0.0d0,0.0d0)
           do 30 i=1,nr-1
              zw=(1.0d0,0.0d0)
              do 25 j=1,nr-1
                 if (j.eq.i) go to 25
                 zw=zw*(z-zo(j))
25            continue
30            zq=zq+zw
           zgd=(zd-zq*zfd)/zp
           z=z-zfd/zgd
           w0=w
           w=cdabs(z)
           if (it.le.50.and.dabs((w-w0)/w).gt.1.0d-12) go to 15
35         zo(nr)=z
        return
        end


        subroutine cfc(z,zf,zd)
c
c       =========================================================
c       purpose: compute complex fresnel integral c(z) and c'(z)
c       input :  z --- argument of c(z)
c       output:  zf --- c(z)
c                zd --- c'(z)
c       =========================================================
c
        implicit double precision (e,p,w)
        implicit complex *16 (c,s,z)
        eps=1.0d-14
        pi=3.141592653589793d0
        w0=cdabs(z)
        zp=0.5d0*pi*z*z
        zp2=zp*zp
        z0=(0.0d0,0.0d0)
        if (z.eq.z0) then
           c=z0
        else if (w0.le.2.5) then
           cr=z
           c=cr
           do 10 k=1,80
              cr=-.5d0*cr*(4.0d0*k-3.0d0)/k/(2.0d0*k-1.0d0)
     &          /(4.0d0*k+1.0d0)*zp2
              c=c+cr
              wa=cdabs(c)
              if (dabs((wa-wa0)/wa).lt.eps.and.k.gt.10) go to 30
10            wa0=wa
        else if (w0.gt.2.5.and.w0.lt.4.5) then
           m=85
           c=z0
           cf1=z0
           cf0=(1.0d-100,0.0d0)
           do 15 k=m,0,-1
              cf=(2.0d0*k+3.0d0)*cf0/zp-cf1
              if (k.eq.int(k/2)*2) c=c+cf
              cf1=cf0
15            cf0=cf
           c=cdsqrt(2.0d0/(pi*zp))*cdsin(zp)/cf*c
        else
           cr=(1.0d0,0.0d0)
           cf=(1.0d0,0.0d0)
           do 20 k=1,20
              cr=-.25d0*cr*(4.0d0*k-1.0d0)*(4.0d0*k-3.0d0)/zp2
20            cf=cf+cr
           cr=1.0d0/(pi*z*z)
           cg=cr
           do 25 k=1,12
              cr=-.25d0*cr*(4.0d0*k+1.0d0)*(4.0d0*k-1.0d0)/zp2
25            cg=cg+cr
           c=.5d0+(cf*cdsin(zp)-cg*cdcos(zp))/(pi*z)
        endif
30      zf=c
        zd=cdcos(0.5*pi*z*z)
        return
        end


        subroutine cfs(z,zf,zd)
c
c       =========================================================
c       purpose: compute complex fresnel integral s(z) and s'(z)
c       input :  z  --- argument of s(z)
c       output:  zf --- s(z)
c                zd --- s'(z)
c       =========================================================
c
        implicit double precision (e,p,w)
        implicit complex *16 (c,s,z)
        eps=1.0d-14
        pi=3.141592653589793d0
        w0=cdabs(z)
        zp=0.5d0*pi*z*z
        zp2=zp*zp
        z0=(0.0d0,0.0d0)
        if (z.eq.z0) then
           s=z0
        else if (w0.le.2.5) then
           s=z*zp/3.0d0
           cr=s
           do 10 k=1,80
              cr=-.5d0*cr*(4.0d0*k-1.0d0)/k/(2.0d0*k+1.0d0)
     &          /(4.0d0*k+3.0d0)*zp2
              s=s+cr
              wb=cdabs(s)
              if (dabs(wb-wb0).lt.eps.and.k.gt.10) go to 30
10            wb0=wb
        else if (w0.gt.2.5.and.w0.lt.4.5) then
           m=85
           s=z0
           cf1=z0
           cf0=(1.0d-100,0.0d0)
           do 15 k=m,0,-1
              cf=(2.0d0*k+3.0d0)*cf0/zp-cf1
              if (k.ne.int(k/2)*2) s=s+cf
              cf1=cf0
15            cf0=cf
           s=cdsqrt(2.0d0/(pi*zp))*cdsin(zp)/cf*s
        else
           cr=(1.0d0,0.0d0)
           cf=(1.0d0,0.0d0)
           do 20 k=1,20
              cr=-.25d0*cr*(4.0d0*k-1.0d0)*(4.0d0*k-3.0d0)/zp2
20            cf=cf+cr
           cr=1.0d0/(pi*z*z)
           cg=cr
           do 25 k=1,12
              cr=-.25d0*cr*(4.0d0*k+1.0d0)*(4.0d0*k-1.0d0)/zp2
25            cg=cg+cr
           s=.5d0-(cf*cdcos(zp)+cg*cdsin(zp))/(pi*z)
        endif
30      zf=s
        zd=cdsin(0.5*pi*z*z)
        return
        end
