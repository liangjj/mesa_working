*deck rmat
      subroutine rmat(rmta,drmta,rmte,drmte,hfbox,eig,rbox,energy,v0,
     1                l,ncon,nen,prnt)
      implicit integer (a-z)
      real*8 hfbox, eig, energy, rmta, drmta, rmte, drmte, v0
      real*8 rbox, beta, dbeta, alpha, dalpha, sn, cn, dsn, dcn
      real*8 psi, dpsi, psidot, dpsidot
      logical prnt
      dimension eig(ncon), hfbox(ncon), energy(nen)
      common /io/ inp, iout
      if (prnt) then
          write(iout,1)
      endif
      do 10 ien=1,nen
         alpha=2.d0*(-v0+energy(ien))
         alpha=sqrt(alpha)
         dalpha=1.d0/alpha
         beta=-2.d0*energy(ien)
         beta=sqrt(beta)
         dbeta=-1.d0/beta
         sn=sin(alpha*rbox)
         cn=cos(alpha*rbox)
         dsn=rbox*cn*dalpha
         dcn=-rbox*sn*dalpha
         if ( l.eq.0 ) then
              psi=sn
              psidot=alpha*cn
              dpsi=dsn
              dpsidot=cn*dalpha+dcn*alpha
         else
             psi=sn/(alpha*rbox) - cn
             psidot=cn/rbox -sn/(alpha*rbox*rbox) + alpha*sn
             dpsi=dsn/(alpha*rbox) - sn*dalpha/(rbox*alpha*alpha) - dcn
             dpsidot=dcn/rbox - dsn/(alpha*rbox*rbox) + 
     1               sn*dalpha/(alpha*alpha*rbox*rbox) +dalpha*sn +
     2               alpha*dsn
         endif 
         rmte=psi/psidot         
         drmte=dsn/psidot -dpsidot*(psi/(psidot*psidot)) 
         rmta=0.d0
         drmta=0.d0
         do 20 i=1,ncon
            rmta=rmta+hfbox(i)*hfbox(i)/(eig(i)-energy(ien))
            drmta=drmta+hfbox(i)*hfbox(i)/
     1                     ((eig(i)-energy(ien))*
     2                     (eig(i)-energy(ien)))         
   20    continue
         rmta=.5d0*rmta
         drmta=.5d0*drmta
         if (prnt) then
             write(iout,2) energy(ien), rmte, rmta, drmte, drmte
         endif    
   10 continue            
    1 format(/,3x,'    energy    ','    ex. rmat   ','    app. rmat    '
     1             ,'   ex. drmat   ','  app. drmat  ')
    2 format(1x,5e15.8) 
      return
      end
