*deck @(#)pdm88.f	5.1  11/6/94
      subroutine pdm88(site,vo,zm,isel,amat,avec,der1,
     $                 zmag,zl,ident,dipole,sdev,xsdev,label,
     $                 nsite,nq,np,nbond,nsave,prnt) 
c***begin prologue     pdm88.f
c***date written       880712
c***revision date      11/6/94
c   12 february, 1994  rlm at lanl
c
c***keywords
c***author             williams, d.e.
c***source             @(#)pdm88.f	5.1   11/6/94
c***purpose
c***description
c
c  variables:
c  ident(i,j)  i=1 has label number of first bond defining atom
c              i=2 has label number of second bond defining atom
c  isel(i,l)   selection integers for zm (0=nosel, 1=sel)
c  nbegin      normally 1, to indicate beginning of atom list or bond list
c                 if greater than 1, indicates the beginning of a bond list
c                 following a saved atom list
c  label(j)    symbol for element (h,c,n,o,f,si,p,s,cl) followed
c                 by sequence number in input
c  site(m,i)   x,y,z (cartesian) for the ith site, m=1,3
c  vo(m,j)     x,y,z,vo for the jth grid point, m=1,4
c                 converted to cartesian during input
c  zl(i,j)     direction cosines of bond j
c  zm(i,l)     lth multipole (variable) of the ith site
c  zmag(i)     magnitude of dipole on atom i
c  amat(jk)    lhs matrix (upper triangle with initial work vector),
c                 converted to inverse and correlation matrix later
c  avec(m)     rhs side of linear equations
c  xsdev(i,l)  esd's of zm
c
c***references
c
c***routines called
c
c***end prologue       pdm88.f
      implicit none
c     --- input variables -----
      integer nsite,np,nq,nbond,nsave
      logical prnt
c     --- input arrays (unmodified) ---
      integer ident(2,nsite)
      integer isel(nsite,10)
      real*8 site(3,nsite)
      real*8 vo(4,np)
      real*8 zl(3,nsite)
      character*5 label(nsite)
c     --- input arrays (scratch) ---
      real*8 avec(nq),der1(nq),sdev(nq),dipole(3)
      real*8 zmag(nsite)
      real*8 amat(nq+nq*(nq+1)/2)
      real*8 xsdev(nsite,10)
c     --- output arrays ---
      real*8 zm(nsite,10)
c     --- output variables ---
c     --- scratch arrays ---
c     --- local variables ---
      integer npar,ip,ivec,i,j,ij,l,l1,isite,ng
      integer k,jk,m,ia,ib,nqq
      integer inp,iout
      real*8 sigma,sigrel,sig2,vc,diff
      real*8 vcalc,vcalcr
      real*8 diptot,debye
      real*8 deriv,clij,zero
c
      parameter (zero=0.0d+00)
c
      common/io/inp,iout
c
 1000 format(/' fits to individual potential points'/
     $'    ip      x         y         z         vc        vo       diff
     $'  )
 1010 format(x,i5,6f10.4)
 1020 format(x/'  sigma='f8.2,'  kj/mol'/' sigrel='f8.2,'  percent')
 1030 format(x/
     $' site         q           ux          uy          uz      magnit
     $ude'/(x,a5,x,4(f10.4,i2),f12.4/
     $     '  esd',2x,    4(f10.4,2x)))
 1040 format(x/
     $' site        qxx         qyy         qzz         qxy         qxz
     $        qyz'/(x,a5,x,6(f10.4,i2)/
     $             '  esd',2x,    6(f10.4,2x)))
 1050 format(x/' bond dipole model, no restriction on direction'/
     $'    bond            ux          uy          uz       magnitude')
 1060 format(x,a5,x,a5,4f12.4/
     $      '  esd',7x,      3f12.4)
 1070 format(' dipole components'3f12.4/' total dipole'f12.4/
     $                                  ' (in debyes) 'f12.4)
 1080 format(x/' bond/atomic dipole model, restricted to direction zl'/
     $'    site            zlx         zly         zlz        value')
 1090 format(x,a5,6x,4f12.4/
     $      '  esd',43x,f12.4)
 1100 format(x,a5,x,a5,4f12.4/
     $      '  esd',43x,       f12.4)
 1110 format(' dipole components'3f12.4/' total dipole'f12.4/
     $                                  ' (in debyes) 'f12.4)
c
c     --- zero matrix and vector
      nqq=nq*(nq+1)/2
      call rzero(amat(nq+1),nqq)
      call rzero(avec,nq)
c
c     --- begin lhs, sum over potential points
      npar=10*nsite
      do 40 ip=1,np
c        --- sum over variables
         ivec=0
         do 10 i=1,npar
            l1=(i-1)/nsite+1
c           --- skip qxx
            if (l1.ne.5) then
               isite=mod(i-1,nsite)+1
c              --- check if this parameter is selected (isel=1)
               if (isel(isite,l1).eq.1) then
                  ivec=ivec+1
c                 --- get derivative
                  if (nbond.lt.2) then
c                    --- normal calculation, nbond=0 or nbond=1
                     der1(ivec)=clij(site(1,isite),vo(1,ip),l1)
                  else
c                    --- restricted calculation, nbond=2
                     der1(ivec)=
     $                  deriv(zl(1,isite),site(1,isite),vo(1,ip))
                  endif
               endif
            endif
   10    continue
c
c        --- finish up matrix, and rhs
         ij=nq
         do 30 i=1,nq
            do 20 j=i,nq
               ij=ij+1
               amat(ij)=amat(ij)+der1(i)*der1(j)
   20       continue
            avec(i)=avec(i)+vo(4,ip)*der1(i)
   30    continue
   40 continue
c
c     --- solve the linear equations, solution comes back in
c         beginning of amat.
      call symlin(amat,nq,avec,ng)
      if (ng.ne.0) then
         call lnkerr('m620: singular matrix')
      endif
c
c     --- distribute variables
      ivec=0
      call rzero(zm,10*nsite)
      do 60 l=1,10
         do 50 isite=1,nsite
            if (isel(isite,l).eq.1) then
               ivec=ivec+1
               zm(isite,l)=amat(ivec)
            endif
   50    continue
   60 continue
c
c     --- get xx quadrupole component
      do 70 isite=1,nsite
         zm(isite,5)=-zm(isite,6)-zm(isite,7)
   70 continue
c
c     --- standard deviations
      sigma=zero
      sigrel=zero
      if (prnt) then
         write(iout,1000)
      endif
      do 80 ip=1,np
         if (nbond.lt.2) then
c           --- normal calculation
            vc=vcalc(nsite,zm,site,vo(1,ip))
         else if (nbond.eq.2) then
c           --- restricted calculation
            vc=vcalcr(nsite,zm,zl,site,vo(1,ip))
         endif
         diff=vc-vo(4,ip)
         if (prnt) then
            write(iout,1010) ip,(vo(i,ip),i=1,3),vc,vo(4,ip),diff
         endif
         sigma=sigma+(vc-vo(4,ip))**2
         sigrel=sigrel+(vo(4,ip))**2
   80 continue
c
      sig2=sigma/(float(np)-float(nq))
      sigma=sqrt(sigma/float(np))
      sigrel=sigrel/float(np)
      sigrel=100.0d0*sigma/(sqrt(sigrel))
      write(iout,1020) sigma,sigrel
c
c     --- convert units from kj/mol to electron charge
c         get esd's
c         scale inverse matrix with sig2
      jk=nq+1
      do 100 j=1,nq
         do 90 k=j,nq
            amat(jk)=amat(jk)*sig2
            if (j.eq.k) sdev(j)=sqrt(amat(jk))
            jk=jk+1
   90    continue
  100 continue
c
c     --- convert to correlation matrix
      jk=nq+1
      do 120 j=1,nq
         do 110 k=j,nq
            amat(jk)=amat(jk)/(sdev(j)*sdev(k))
            jk=jk+1
  110    continue
  120 continue
c
c     --- scale esd's
      do 130 j=1,nq
         sdev(j)=sdev(j)/1389.3655d0
  130 continue
c     --- length unit is the angstrom
      do 150 isite=1,nsite
         do 140 l=1,10
            zm(isite,l)=zm(isite,l)/1389.3655d0
  140    continue
         zmag(isite)=sqrt(zm(isite,2)**2
     $              +zm(isite,3)**2
     $              +zm(isite,4)**2)
  150 continue
c
c     --- distribute esd's
      ivec=0
      call rzero(xsdev,10*nsite)
      do 170 l=1,10
         do 160 isite=1,nsite
            if (isel(isite,l).eq.1) then
               ivec=ivec+1
               xsdev(isite,l)=sdev(ivec)
            endif
  160    continue
  170 continue
c
c     --- print the results
      if (nbond.eq.0) then
         write(iout,1030) 
     $      (label(isite),(zm(isite,l),isel(isite,l),l=1,4),
     $      zmag(isite),(xsdev(isite,m),m=1,4),isite=1,nsite)
         write(iout,1040) 
     $      (label(isite),(zm(isite,l),isel(isite,l),l=5,10),
     $      (xsdev(isite,l),l=5,10),isite=1,nsite)
         if(prnt) call outrho(amat,nq)
      else if (nbond.eq.1) then
         write(iout,1050)
         call rzero(dipole,3)
         do 200 isite=1,nsite
            ia=ident(1,isite)
            ib=ident(2,isite)
            do 190 i=1,3
               dipole(i)=dipole(i)+zm(isite,i+1)
  190       continue
            write(iout,1060) label(ia),label(ib),
     $                       (zm(isite,l),l=2,4),zmag(isite),
     $                       (xsdev(isite,l),l=2,4)
  200    continue
         diptot=sqrt(dipole(1)**2+dipole(2)**2+dipole(3)**2)
         debye=diptot*4.803241d0
         write(iout,1070) (dipole(i),i=1,3),diptot,debye
         if(prnt) call outrho(amat,nq)
      else if (nbond.eq.2) then
         write(iout,1080)
         if (nsave.ne.0) then
            call rzero(dipole,3)
c           --- output saved atomic sites
            do 220 isite=1,nsave
               write(iout,1090) label(isite),(zl(i,isite),i=1,3),
     $                          zm(isite,1),xsdev(isite,1)
  220       continue
         endif
c        --- output bond sites
         do 230 isite=nsave+1,nsite
            ia=ident(1,isite)
            ib=ident(2,isite)
            write(iout,1100) label(ia),label(ib),
     $                       (zl(l,isite),l=1,3),zm(isite,1),
     $                       xsdev(isite,1)
  230    continue
c        --- accumulate resultant dipole moment, all sites
         do 250 isite=1,nsite
            do 240 i=1,3
               dipole(i)=dipole(i)+zl(i,isite)*zm(isite,1)
  240       continue
  250    continue
c
c
         diptot=sqrt(dipole(1)**2+dipole(2)**2+dipole(3)**2)
         debye=diptot*4.803241d0
         write(iout,1110) (dipole(i),i=1,3),diptot,debye
         if(prnt) call outrho(amat,nq)
      endif
c
c
      return
      end
