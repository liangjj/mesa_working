c $Header: shells.f,v 1.1 92/12/31 15:23:30 bis Exp $
*deck @(#)shells.f	1.1 9/9/91
      subroutine shells(rpt,thpt,phpt,wtr,wtth,wtph,ns,nrmax,
     1                  nthmax,nphmax,which,prnt)
      parameter (numshl=50, nptdim=100)
      implicit real*8 (a-h,o-z) 
      real*8 nmulphi,nmulth
      logical prnt
      character*20 itype, jtype, ktype
      character*2 itoc, chrv
      character*(*) which
      dimension rpt(nrmax,ns), thpt(nthmax,ns), phpt(nphmax,ns)
      dimension wtr(nrmax,ns), wtth(nthmax,ns), wtph(nphmax,ns)
      common/spheri/ nshell,nr(numshl),nthell(numshl), nphell(numshl),
     1               ntheta(numshl,numshl), nphi(numshl,numshl)
      common/spherr/ r(numshl),theta(numshl,numshl),phi(numshl,numshl)
      common /gauss/ x(nptdim,3),w(nptdim,3),scr(nptdim),dummy(2)
      common/io/inp,iout
      common/quadtp/ itype, jtype, ktype
      pi=3.14159265358d+00
      fac=pi/180.d+00
      do 10 i=1,numshl
         do 20 j=1,numshl
            phi(i,j)=fac*phi(i,j)
            theta(i,j)=theta(i,j)*fac
   20    continue
   10 continue
      fac1=theta(nthell(1)+1,1)-theta(1,1)
      fac2=phi(nphell(1)+1,1)-phi(1,1)
      do 30 i=2,nshell
         fac3=theta(nthell(i)+1,i)-theta(1,i)
         fac4=phi(nphell(i)+1,i)-phi(1,i)
         if(fac1.ne.fac3.or.fac2.ne.fac4) then
            call lnkerr(' theta/phi boundary inconsistent')
         endif
   30 continue
c     multiplicative factors for theta and phi integration if symmetry
c     allows integrals to be calculated on subdomains.
      nmulphi=2.d+00*pi/fac2
      nmulth=pi/fac1
      write(iout,100) nmulth, nmulphi
      icount=0
c
c open a loop on shells
c
      do 40 i=1,nshell
c     get radial quadrature points and weights for this shell
         ampr=(r(i+1)-r(i))/2.d+00
         abpr=(r(i+1)+r(i))/2.d+00
         if (i.eq.1) then
             call gaussq(itype,nr(i),0.d+00,0.d+00,0,dummy,scr,
     1                   x(1,1),w(1,1))
         elseif (nr(i).ne.nr(i-1)) then
             call gaussq(itype,nr(i),0.d+00,0.d+00,0,dummy,scr,
     1                   x(1,1),w(1,1))
         endif
         do 50 j=1,nr(i)
            numr=numr+1
            rpt(j,i)=ampr*x(j,1)+abpr
            wtr(j,i)=rpt(j,i)*rpt(j,i)*ampr*w(j,1)
   50    continue
         numthe=0
         do 60 j=1,nthell(i)
            ampt=(cos(theta(j,i))-cos(theta(j+1,i)))/2.d+00
            abpt=(cos(theta(j,i))+cos(theta(j+1,i)))/2.d+00
            if (j.eq.1) then
                call gaussq(jtype,ntheta(j,i),0.d+00,0.d+00,0,
     1                      dummy,scr,x(1,2),w(1,2))
            elseif (ntheta(j,i).ne.ntheta(j-1,i)) then
                call gaussq(jtype,ntheta(j,i),0.d+00,0.d+00,0,
     1                      dummy,scr,x(1,2),w(1,2))
            endif
            do 70 k=1,ntheta(j,i)
               numthe=numthe+1
               thpt(numthe,i)=ampt*x(k,2)+abpt
               wtth(numthe,i)=nmulth*ampt*w(k,2)
   70       continue
   60    continue
         numphi=0
         do 80 j=1,nphell(i)
            ampp=(phi(j+1,i)-phi(j,i))/2.d+00
            abpp=(phi(j+1,i)+phi(j,i))/2.d+00
            if(j.eq.1) then
               call gaussq(ktype,nphi(j,i),0.d+00,0.d+00,0,dummy,
     1                     scr,x(1,3),w(1,3))
            elseif (nphi(j,i).ne.nphi(j-1,i)) then
               call gaussq(ktype,nphi(j,i),0.d+00,0.d+00,0,dummy,
     1                     scr,x(1,3),w(1,3))
            endif
            do 90 k=1,nphi(j,i)
               numphi=numphi+1
               phpt(numphi,i)=ampp*x(k,3)+abpp
               wtph(numphi,i)=nmulphi*ampp*w(k,3)
   90       continue
   80    continue
         chrv=itoc(i)
         if (which.eq.'all') then
             call iosys ('write real "radial points shell-'//chrv//
     1                   '" to lamdat',nr(i),rpt(1,i),0,' ')
             call iosys ('write real "radial weights shell-'//chrv//
     1                   '" to lamdat',nr(i),wtr(1,i),0,' ')
             call iosys ('write real "theta points shell-'//chrv//
     1                   '" to lamdat',numthe,thpt(1,i),0,' ')
             call iosys ('write real "theta weights shell-'//chrv//
     1                   '" to lamdat',numthe,wtth(1,i),0,' ')
             call iosys ('write real "phi points shell-'//chrv//
     1                   '" to lamdat',numphi,phpt(1,i),0,' ')
             call iosys ('write real "phi weights shell-'//chrv//
     1                   '" to lamdat',numphi,wtph(1,i),0,' ')
         else
             call iosys ('write real "radial points for atom-'//which//
     1                   '-shell-'//chrv//'" to lamdat',nr(i),
     2                    rpt(1,i),0,' ')
             call iosys ('write real "radial weights for atom-'//which//
     1                   '-shell-'//chrv//'" to lamdat',nr(i),
     2                    wtr(1,i),0,' ')
             call iosys ('write real "theta points for atom-'//which//
     1                   '-shell-'//chrv//'" to lamdat',numthe,
     2                   thpt(1,i),0,' ')
             call iosys ('write real "theta weights for atom-'//which//
     1                   '-shell-'//chrv//'" to lamdat',numthe,
     2                    wtth(1,i),0,' ')
             call iosys ('write real "phi points for atom-'//which//
     1                    '-shell-'//chrv//'" to lamdat',numphi,
     2                     phpt(1,i),0,' ')
             call iosys ('write real "phi weights for atom-'//which//
     1                   '-shell-'//chrv//'" to lamdat',numphi,
     2                    wtph(1,i),0,' ')
         endif
         if (prnt) then
             write(iout,*) '     data on points for shell = ',i
             write(iout,*)
             write(iout,*) 'radial points and weights'
             write(iout,110) (rpt(ii,i),ii=1,nr(i))
             write(iout,110) (wtr(ii,i),ii=1,nr(i))
             write(iout,*)
             write(iout,*) 'cos(theta) points and weights'
             write(iout,110) (thpt(ii,i),ii=1,numthe)
             write(iout,110) (wtth(ii,i),ii=1,numthe)
             write(iout,*)
             write(iout,*) 'phi points and weights'
             write(iout,110) (phpt(ii,i),ii=1,numphi)
             write(iout,110) (wtph(ii,i),ii=1,numphi)
         endif
   40 continue
      return
  100 format(' multiplicative factors for theta and phi  :',2f10.5)
  110 format(/4(2x,e15.8))
      end
