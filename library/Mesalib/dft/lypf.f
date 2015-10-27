*deck @(#)lypf.f	5.2  2/5/95
      subroutine lypf(ngrid,mxgrd,dengrida,dengridb,ga,gb,gab,omega,
     $     delta,t,derivs,a,b,c,d,fout)
c***begin prologue     lypf.f
c***date written       930510     (yymmdd)
c***revisiondate       2/5/95
c
c***keywords           xm602, link 602, DFT, LYP, gradient corrected
c
c***author             RUSSO, thomas v.    (lanl)
c***source             @(#)lypf.f	5.2   2/5/95
c***description        computes the value of the LYP functional
c                      given the value of the density on a grid of points.
c***description2
c                      THIS VERSION USES TOO MUCH DADGUM TEMP SPACE!
c                      C'mon, get with the program, reduce, reuse, recycle!
c
c***references
c        Johnson, B. G. et al., J. Chem. Phys 98(7),5612
c***routines called
c                      
c
c***end prologue       lypf.f
c
c    ----------
c death to the bloated lackey of imperialist FORTRASH, implicit typing,
c the bane of mine existence and my mortal enemy
c
      implicit none
c
c --- input variables ---
c
c derivs>0 do both functional and derivs
c derivs=0 do only functional
c derivs<0 do only derivs
      integer ngrid,derivs,mxgrd
      real*8 a,b,c,d
c
c --- input arrays (unmolested) ---
c
      real*8 dengrida(ngrid),dengridb(ngrid),ga(ngrid),gb(ngrid),
     $     gab(ngrid)
c
c --- input arrays (scribbled all over then discarded) ---
c
      real*8 omega(ngrid),delta(ngrid),t(mxgrd,*)
c --- input arrays (modified)
c --- Rule Psix ---
c (there is No Rule Psix)
c
c --- output arrays ---
      real*8 fout(mxgrd,*)
c
c
c --- sample the local cuisine ---      
c
      integer i
      real*8 twoelev3,pif,one,two,three,four,seven,nine,eleven,forty7,
     $     un3,eight3
      real*8 pi
      integer inp,iout
      common /io/inp,iout
c
      one=1.0d0
      two=2.0d0
      three=3.0d0
      four=4.0d0
      seven=7.0d0
      nine=9.0d0
      eleven=11.0d0
      forty7=47.0d0
      pi=four*atan(one)
      twoelev3=two**(eleven/three)*(three/10.0d0)
      pif=(three*pi*pi)**(two/three)
      un3=one/three
      eight3=8.0d0/three
c
      call rzero(t,mxgrd*10)
      call rzero(delta,mxgrd)
      call rzero(omega,mxgrd)
      do 10 i=1,ngrid
         t(i,1)=dengrida(i)+dengridb(i)
 10   continue 
c
c
      do 15 i=1,ngrid
            t(i,2)=t(i,1)**(-un3)
            t(i,3)=one/(one+d*t(i,2))
            delta(i)=(c+d*t(i,3))*t(i,2)
            omega(i)=exp(-c*t(i,2))*t(i,3)*(t(i,2)**eleven)
            t(i,4)=dengrida(i)**eight3+dengridb(i)**eight3
            t(i,5)=dengrida(i)*dengridb(i)
 15   continue 
c
c dLYP/dgamma, but don't stick in fout unless derivs nonzero
c
      do 20 i=1,ngrid
            t(i,6)=-a*b*omega(i)*(t(i,5)/nine*(one-three*delta(i)-
     $           (delta(i)-eleven)*dengrida(i)/t(i,1))-
     $           dengridb(i)**two)
            t(i,7)=-a*b*omega(i)*(t(i,5)/nine*(one-three*delta(i)-
     $           (delta(i)-eleven)*dengridb(i)/t(i,1))-
     $           dengrida(i)**two)
            t(i,8)=-a*b*omega(i)*(t(i,5)/nine*(forty7-seven*delta(i))-
     $           four/three*t(i,1)*t(i,1))
 20   continue 
      if (derivs .ne.0)then
         call vadd(fout(1,3),fout(1,3),t(1,6),ngrid)
         call vadd(fout(1,5),fout(1,5),t(1,7),ngrid)
         call vmove(fout(1,6),t(1,8),ngrid)
      endif
c
      if (derivs .ge. 0) then
         do 30 i=1,ngrid
               fout(i,1)=fout(i,1)-four*a*t(i,3)*t(i,5)/t(i,1)
               fout(i,1)=fout(i,1)-twoelev3*pif*a*b*omega(i)*
     $              t(i,5)*t(i,4)
               fout(i,1)=fout(i,1)+t(i,6)*ga(i)+t(i,7)*gb(i)+
     $              t(i,8)*gab(i)
 30      continue 
      endif
c
c omega' and delta'
c
      if (derivs .ne. 0) then
         do 40 i=1,ngrid
               t(i,9)=-t(i,2)**four*omega(i)/three*(eleven/t(i,2)-c-
     $              d*t(i,3))
               t(i,10)=(d*d*t(i,2)**5*t(i,3)**2-delta(i)/t(i,1))/three
 40      continue 
c
c do the stuff that doesn't need the second partials
c
         do 50 i=1,ngrid
               if (dengrida(i).ne.0)
     $          fout(i,2)=fout(i,2)-4*a*t(i,3)*t(i,5)/t(i,1)*
     $              (d*t(i,2)**4*t(i,3)/three+one/dengrida(i)-
     $              one/t(i,1))
               fout(i,2)=fout(i,2)-twoelev3*pif*a*b*(t(i,9)*t(i,5)*
     $              t(i,4)+omega(i)*dengridb(i)*
     $              (eight3*dengrida(i)**(eight3)+t(i,4)))
               if (dengridb(i).ne.0)
     $          fout(i,4)=fout(i,4)-4*a*t(i,3)*t(i,5)/t(i,1)*
     $              (d*t(i,2)**4*t(i,3)/three+one/dengridb(i)-
     $              one/t(i,1))
               fout(i,4)=fout(i,4)-twoelev3*pif*a*b*(t(i,9)*t(i,5)*
     $              t(i,4)+omega(i)*dengrida(i)*
     $              (eight3*dengridb(i)**(eight3)+t(i,4)))
 50      continue 
c
c now calculate the mixed second partials of LYP and finish the deed.
c
         do 55 i=1,ngrid
            if (omega(i) .ne.0) then
               fout(i,2)=fout(i,2)+t(i,9)/omega(i)*(ga(i)*t(i,6)+
     $              gb(i)*t(i,7)+gab(i)*t(i,8))
               fout(i,4)=fout(i,4)+t(i,9)/omega(i)*(ga(i)*t(i,6)+
     $              gb(i)*t(i,7)+gab(i)*t(i,8))
            endif
 55      continue 
         do 60 i=1,ngrid
c
c the dLYP/d(rho)d(gaa) terms
c
            fout(i,2)=fout(i,2)+ga(i)*(
     $           -a*b*omega(i)*(dengridb(i)/nine*
     $           (one-three*delta(i)-(delta(i)-eleven)*dengrida(i)
     $           /t(i,1))-
     $           t(i,5)/nine*((three+dengrida(i)/t(i,1))*t(i,10)+
     $           (delta(i)-eleven)*dengridb(i)/t(i,1)**2)))
            fout(i,4)=fout(i,4)+ga(i)*(
     $           -a*b*omega(i)*(dengrida(i)/nine*
     $           (one-three*delta(i)-(delta(i)-eleven)*dengrida(i)
     $           /t(i,1))-
     $           t(i,5)/nine*((three+dengrida(i)/t(i,1))*t(i,10)-
     $           (delta(i)-eleven)*dengrida(i)/t(i,1)**2)-
     $           two*dengridb(i)))
 60      continue 
         do 70 i=1,ngrid
c
c the dLYP/d(rho)d(gbb) terms
c
            fout(i,2)=fout(i,2)+gb(i)*(
     $           -a*b*omega(i)*(dengridb(i)/nine*
     $           (one-three*delta(i)-(delta(i)-eleven)*dengridb(i)
     $              /t(i,1))-
     $           t(i,5)/nine*((three+dengridb(i)/t(i,1))*t(i,10)-
     $           (delta(i)-eleven)*dengridb(i)/t(i,1)**2)-
     $              2*dengrida(i)))
            fout(i,4)=fout(i,4)+gb(i)*(
     $           -a*b*omega(i)*(dengrida(i)/nine*
     $           (one-three*delta(i)-(delta(i)-eleven)*dengridb(i)
     $           /t(i,1))-
     $           t(i,5)/nine*((three+dengridb(i)/t(i,1))*t(i,10)+
     $           (delta(i)-eleven)*dengrida(i)/t(i,1)**2)))
 70      continue 
c
c now the dLYP/d(rho)d(gab) terms
c
         do 80 i=1,ngrid
            fout(i,2)=fout(i,2)+gab(i)*(
     $           -a*b*omega(i)*(dengridb(i)/nine*
     $           (forty7-seven*delta(i))-seven*t(i,5)/nine*t(i,10)-
     $              eight3*t(i,1)))
            fout(i,4)=fout(i,4)+gab(i)*(
     $           -a*b*omega(i)*(dengrida(i)/nine*
     $           (forty7-seven*delta(i))-seven*t(i,5)/nine*t(i,10)-
     $           eight3*t(i,1)))
 80      continue 
      endif

      return
      end
