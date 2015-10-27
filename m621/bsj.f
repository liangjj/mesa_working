*deck @(#)bsj.f	5.1 11/6/94
      subroutine bsj(x,y,z,q,qd,qq,
     $         xcnt,ycnt,zcnt,rad,
     $         xobs,yobs,zobs,epsilon,npoints,nsrcs,
     $         ncnt,nobs,prnt,aben,deltag2,pomf2,self,
     $         nfocmx,mark,nark,
     $         phi,ex,ey,ez,efgxx,efgxy,
     $         efgxz,efgyy,efgyz,efgzz,
     $         xq,yq,zq,indx,u,
     $         a2,b2,nfocus,xyzch,charge,ifocus,
     $         icindx,ipindx,index1,index2,
     $         index3,toang,toev,tokcal)
c***begin prologue     bsj.f
c***date written       yymmdd
c***revision date      11/6/94
c   august 5, 1994     gjt at lanl
c      adding capability to handle through l=2 in the atom-centered multipole
c      expansion of the potential.
c***keywords
c***author             tawa,greg(lanl)
c***source             @(#)bsj.f	5.1   11/6/94
c***purpose
c***description
c
c
c
c***references
c
c***routines called
c
c***end prologue       bsj.f
      implicit none
c     --- input variables -----
      integer npoints,nsrcs,ncnt,nobs,nfocmx
      logical prnt
      real*8 epsilon
      real*8 aben
      real*8 toang,toev,tokcal
c     --- input arrays (unmodified) ---
      real*8 x(nsrcs),y(nsrcs),z(nsrcs),q(nsrcs)
      real*8 qd(nsrcs,3),qq(nsrcs,6)
      real*8 xcnt(ncnt),ycnt(ncnt),zcnt(ncnt),rad(ncnt)
      real*8 xobs(nobs),yobs(nobs),zobs(nobs)
c     --- input arrays (scratch) ---
      integer ifocus(ncnt)
      integer icindx(npoints),ipindx(npoints)
      integer index2(npoints),index3(nfocmx)
      integer indx(npoints)
      integer mark(ncnt+1),nark(ncnt+1)
      real*8 phi(nobs)
      real*8 ex(nobs),ey(nobs),ez(nobs),efgxx(nobs),efgxy(nobs)
      real*8 efgxz(nobs),efgyy(nobs),efgyz(nobs),efgzz(nobs)
      real*8 u(2*ncnt)
      real*8 a2(nfocmx,nfocmx),b2(nfocmx)
      real*8 xq(npoints),yq(npoints),zq(npoints)
      real*8 rc,r,r2,r3,r5,r7,dot1,dot2x,dot2y,dot2z,dotd,f,f2
      real*8 rqxx,rqxy,rqxz,rqyy,rqyz,rqzz
      real*8 rqx,rqy,rqz
      real*8 dx22,dy22,dz22,dotm
c     --- output arrays ---
      integer index1(nfocmx)
      real*8 xyzch(3,nfocmx),charge(nfocmx)
c     --- output variables ---
      integer nfocus
      real*8 deltag2,pomf2,self
      real*8 e2
c     --- scratch arrays ---
c     --- local variables ---
      integer inp,iout
      integer i,j,kount
      integer k
      integer isrc,iobs
      real*8 zero,one,two,four
      real*8 pi,twopi,fourpi
      real*8 d
      real*8 factor1,dx,dy,dz,dx2,dy2,dz2,d2
      real*8 dot,e,factor,twosigma
      real*8 hudge
      real*8 total
      parameter (hudge=1.0d+10)
      logical timing
      logical debug
c
      data zero/0.0d+00/,one/1.0d+00/,two/2.0d+00/,four/4.0d+00/
      save zero,one,two,four
c
      parameter (timing=.false.)
      parameter (debug=.false.)
c     parameter (debug=.true.)
c
      common/io/inp,iout
c
 1000 format(5x,'solvent free energy    ',e15.9)
 1010 format(5x,'potential of mean force',e15.9)
 1020 format(5x,'gas phase energy       ',e15.9)
 1030 format(5x,'solvent self energy    ',e15.9)
c
c     --- set up some constants and times ---
      pi=four*atan(one)
      twopi=two*pi
      fourpi=four*pi
c
c     --- get the surface and focus points.
      call surface(xyzch,xcnt,ycnt,zcnt,rad,
     $             npoints,ncnt,prnt,nfocmx,mark,nark,
     $             xq,yq,zq,indx,u,kount,
     $             nfocus,ifocus,icindx,ipindx,index1,index2,index3)
c
c     --- set up the matrix a2
c         Left transform
      e2 = one
      do 90 j = 1,nfocus
         do 95 k = 1,nfocus
            a2(j,k)=zero
c           area of plac; focus point k
            factor1= (fourpi*rad(icindx(index1(k)))**2)/
     &               float(ipindx(index1(k)))
            do 100 i = 1,kount
               if (index2(i) .eq. j) then
                  if (i .eq.index1(k)) then
                     a2(j,k) = a2(j,k) +
     &                        one +
     &                        (one-one/sqrt(float(ipindx(i))))*
     &                        (e2-epsilon)/(two*epsilon)
                  else
                     dx=xq(i)-xcnt(icindx(i))
                     dy=yq(i)-ycnt(icindx(i))
                     dz=zq(i)-zcnt(icindx(i))
                     d = sqrt(dx**2+dy**2+dz**2)
                     dx2 = xq(i)-xq(index1(k))
                     dy2 = yq(i)-yq(index1(k))
                     dz2 = zq(i)-zq(index1(k))
                     d2 = sqrt(dx2*dx2+dy2*dy2+dz2*dz2)
                     dot = dx*dx2+dy*dy2+dz*dz2
                     a2(j,k) = a2(j,k)
     $                       + ((epsilon -e2)/(fourpi*epsilon))
     &                       * dot*factor1/(d*d2*d2*d2)
                  endif
               endif
  100       continue
  95     continue
  90  continue
c
c     --- print out matrix a2
      if (debug) then
         write(iout,'(''Matrix for normal E             '')')
         call matout(a2,nfocmx,nfocmx,nfocus,nfocus,iout)
      endif
c
c     --- now compute the rhs vector
c         first convert the dipole and quadrupole components to atomic units
c
      do 123 i=1,nsrcs
         qd(i,1)=qd(i,1)/toang
         qd(i,2)=qd(i,2)/toang
         qd(i,3)=qd(i,3)/toang
         qq(i,1)=qq(i,1)/(toang*toang)
         qq(i,2)=qq(i,2)/(toang*toang)
         qq(i,3)=qq(i,3)/(toang*toang)
         qq(i,4)=qq(i,4)/(toang*toang)
         qq(i,5)=qq(i,5)/(toang*toang)
         qq(i,6)=qq(i,6)/(toang*toang)
  123 continue
      if (debug) then
         write(iout,'(''Right hand side vector for normal E'')')
      endif
      do 105 j = 1,nfocus
         b2(j) = zero
         do 110 i = 1,kount
            if (index2(i) .ne. 0 .and. index2(i) .eq. j) then
               dx = xq(i) -xcnt(icindx(i))
               dy = yq(i) -ycnt(icindx(i))
               dz = zq(i) -zcnt(icindx(i))
               rc = sqrt(dx**2+dy**2+dz**2)
               do 115 isrc = 1,nsrcs
                  dx2 = xq(i)-x(isrc)
                  dy2 = yq(i)-y(isrc)
                  dz2 = zq(i)-z(isrc)
                  dx22=dx2*dx2
                  dy22=dy2*dy2
                  dz22=dz2*dz2
                  r = sqrt(dx2**2+dy2**2+dz2**2)
                  r2=r*r
                  r3=r*r*r
                  r5=r*r*r*r*r
                  r7=r*r*r*r*r*r*r
                  dotm = dx*dx2+dy*dy2+dz*dz2
                  dot1=qd(isrc,1)*dx2+qd(isrc,2)*dy2+qd(isrc,3)*dz2
                  dot2x=3.0d0*dot1*dx2/r5-qd(isrc,1)/r3
                  dot2y=3.0d0*dot1*dy2/r5-qd(isrc,2)/r3
                  dot2z=3.0d0*dot1*dz2/r5-qd(isrc,3)/r3
                  dotd=dx*dot2x+dy*dot2y+dz*dot2z
                  f = one
                  f2 = one
c                 --- x component, quadrupole component of electric field
                  rqxx=(qq(isrc,1)*dx2/(two*r5))*
     &                 (5.d0*dx22/r2-two*f2)
                  rqyy=5.d0*qq(isrc,2)*dy22*dx2/(two*r7)
                  rqzz=5.d0*qq(isrc,3)*dz22*dx2/(two*r7)
                  rqxy=(3.d0*qq(isrc,4)*dy2/r5)*f*
     &                 (5.d0*dx22/r2-one*f2)
                  rqxz=(3.d0*qq(isrc,5)*dz2/r5)*f*
     &                 (5.d0*dx22/r2-one*f2)
                  rqyz=f*15.d0*qq(isrc,6)*dy2*dz2*dx2/r7
                  rqx=(dx/rc)*(rqxx+rqyy+rqzz+rqxy+rqxz+rqyz)
c                 --- y component, quadrupole component of electric field
                  rqxx=5.d0*qq(isrc,1)*dx22*dy2/(two*r7)
                  rqyy=(qq(isrc,2)*dy2/(two*r5))*
     &                 (5.d0*dy22/r2-two*f2)
                  rqzz=5.d0*qq(isrc,3)*dz22*dy2/(two*r7)
                  rqxy=(3.d0*qq(isrc,4)*dx2/r5)*f*
     &                 (5.d0*dy22/r2-1.d0*f2)
                  rqxz=f*15.d0*qq(isrc,5)*dy2*dz2*dx2/r7
                  rqyz=(3.d0*qq(isrc,6)*dz2/r5)*f*
     &                 (5.d0*dy22/r2-one*f2)
                  rqy=(dy/rc)*(rqxx+rqyy+rqzz+rqxy+rqxz+rqyz)
c                 --- z component, quadrupole component of electric field
                  rqxx=5.d0*qq(isrc,1)*dx22*dz2/(two*r7)
                  rqyy=5.d0*qq(isrc,2)*dy22*dz2/(two*r7)
c
c                 BUG FIX 7-7-94 dy2 on next line changed to dz2
                  rqzz=(qq(isrc,3)*dz2/(two*r5))*
     &                 (5.d0*dz22/r2-two*f2)
                  rqxy=f*15.d0*qq(isrc,4)*dy2*dz2*dx2/r7
                  rqxz=(3.d0*qq(isrc,5)*dx2/r5)*f*
     &                 (5.d0*dz22/r2-one*f2)
                  rqyz=(3.d0*qq(isrc,6)*dy2/r5)*f*
     &                 (5.d0*dz22/r2-one*f2)
                  rqz=(dz/rc)*(rqxx+rqyy+rqzz+rqxy+rqxz+rqyz)
c
c                 ---  full right hand side vector = electric field due to solute
c                      charge distribution
                  b2(j)=b2(j)+
     &                  q(isrc)*dotm/(rc*r3)+
     &                  dotd/rc+rqx+rqy+rqz
  115          continue
            endif
  110    continue
         if (debug) then
           write(iout,'(e15.9)') b2(j)
         endif
         b2(j) = b2(j)/e2
  105 continue
c
c     --- done with matrix and r-side vector -- now solve
      if(debug) then
         write(iout, 530)
  530    format(/4x,'BOTH MATRIX AND VECTOR ARE NOW READY -- SOLVE!'/)
      endif
      call ludcmp(a2,nfocus,nfocmx,indx,e)
      call lubksb(a2,nfocus,nfocmx,indx,b2)
c
c     --- now calculate the surface charges and the induced potentials
c         at the observation points due to the surface charges
      total = zero
      if(debug) then
         write(iout,*) 'suface charges and positions'
      endif
      do 120 i = 1,nfocus
         factor = (rad(icindx(index1(i)))**2)*(e2-epsilon)/
     &            (epsilon*float(ipindx(index1(i))))
         charge(i)=b2(i)*factor
         total = total + charge(i)
         if(debug) then
            write(iout,
     $              '(f10.5,1x,f10.5,1x,f10.5,1x,e15.9,1x,e15.9,
     &              i10)')
     &              xq(index1(i)),yq(index1(i)),zq(index1(i)),
     &              charge(i),total,index3(i)
         endif
  120 continue
c
c     --- write out the normal component of the electric field vector
c         to be used as an initial guess in an importance sampling
c         monte carlo calculation
c     do 125 i = 1,nfocus
c        write(25,'(f10.5,1x,f10.5,1x,f10.5,1x,e15.9)')
c    &             xq(index1(i)),yq(index1(i)),zq(index1(i)),
c    &             b2(i)
c 125 continue
c     --- calculate the induced potential at observation points due
c         to the point charges on the molecular surface
      call rzero(phi,nobs)
      call rzero(ex,nobs)
      call rzero(ey,nobs)
      call rzero(ez,nobs)
      call rzero(efgxx,nobs)
      call rzero(efgxy,nobs)
      call rzero(efgxz,nobs)
      call rzero(efgyy,nobs)
      call rzero(efgyz,nobs)
      call rzero(efgzz,nobs)
c
      do 130 iobs = 1,nobs
         do 135 i = 1,nfocus
            dx = x(iobs)-xq(index1(i))
            dy = y(iobs)-yq(index1(i))
            dz = z(iobs)-zq(index1(i))
            r  = sqrt(dx*dx+dy*dy+dz*dz)
            r3 = r*r*r
            r5 = r3*r*r
c           --- potential
            phi(iobs) = phi(iobs) + charge(i)/r
c           --- electrie field components
            ex(iobs)  = ex(iobs)  + charge(i)*dx/r3
            ey(iobs)  = ey(iobs)  + charge(i)*dy/r3
            ez(iobs)  = ez(iobs)  + charge(i)*dz/r3
c           --- electrie field gradient components
            efgxx(iobs)=efgxx(iobs)+charge(i)*(1.d0/r3-3.d0*dx*dx/r5)
            efgyy(iobs)=efgyy(iobs)+charge(i)*(1.d0/r3-3.d0*dy*dy/r5)
            efgzz(iobs)=efgzz(iobs)+charge(i)*(1.d0/r3-3.d0*dz*dz/r5)
            efgxy(iobs)=efgxy(iobs)-3.0d0*charge(i)*dx*dy/r5
            efgxz(iobs)=efgxz(iobs)-3.0d0*charge(i)*dx*dz/r5
            efgyz(iobs)=efgyz(iobs)-3.0d0*charge(i)*dy*dz/r5
  135    continue
  130 continue
c
c     ---  calculate the interaction energy of the surface charges
c     The factor of .5 is present in the equation for self in order that
c     interactions are not double counted.
      self = zero
      do 140 i=1,nfocus
         do 145 j=1,nfocus
            if (j .ne. i) then
               dx=xq(index1(i))-xq(index1(j))
               dy=yq(index1(i))-yq(index1(j))
               dz=zq(index1(i))-zq(index1(j))
               d=dx*dx+dy*dy+dz*dz
               self = self + charge(i)*charge(j)/sqrt(d)
            endif
  145    continue
  140 continue
      self = self*0.5d0
c
c     --- calculate the electrostatic component of the solvation
c         free energy and the potential of mean force
      deltag2=zero
      do 155 isrc = 1,nsrcs
         if(debug) then
            write(iout,'(/4x,''Observation points'',t30,3f10.5,1x)')
     &                   xobs(isrc),yobs(isrc),zobs(isrc)
            write(iout,'(4x,''Source points'',t30,3f10.5,1x)')
     &                   x(isrc),y(isrc),z(isrc)
            write(iout,'(4x,''Charges and potentials'',t30,3f10.5,1x)')
     &                   q(isrc),phi(isrc)
         endif
         deltag2 = deltag2
     &            +q(isrc)*phi(isrc)
     &            -(ex(isrc)*qd(isrc,1)
     &            +ey(isrc)*qd(isrc,2)+ez(isrc)*qd(isrc,3))
     &            -(one/6.0d0)*(qq(isrc,1)*efgxx(isrc)
     &            +qq(isrc,2)*efgyy(isrc)
     &            +qq(isrc,3)*efgzz(isrc)
     &            +6.0d0*qq(isrc,4)*efgxy(isrc)
     &            +6.0d0*qq(isrc,5)*efgxz(isrc)
     &            +6.0d0*qq(isrc,6)*efgyz(isrc))
  155 continue
c
c     --- note the "mysterious" factor of 0.5 which 
c     reflects the fact that the interaction is induced.
c
      deltag2=deltag2/two
c     --- calculate the potential of Mean force
      pomf2 = deltag2 + aben
      if(debug) then
         write(iout,540)
  540    format(/2x,'FINAL RESULTS FOR THE OBSERVATION POINTS')
      endif
c
c     note that atomic units have been used throughout. convert to
c     more comfortable ones now.
c     toev(~27.21 eV) converts the computed potential to VOLTS
c     toang(~.52917A) converts bohr to angstroms.
c     tokcal(~627.1) converts hartree to kcal/mole.
c     electric fields are then VOLTS/BOHR, etc...
      if(debug) then
         do 160 iobs = 1, nobs
            write(iout, 550)iobs, xobs(iobs), yobs(iobs), zobs(iobs)
  550	    format(/2x,'observation point',i3
     &         ,' at xobs = ',f6.3,2x,'yobs = ',f6.3,2x,'zobs = ',f6.3)
	    twosigma = zero
	    write(iout,560)	toev*phi(iobs), twosigma
     	    write(iout,570)	(toev/toang)*Ex(iobs)
     &,				(toev/toang)*Ey(iobs)
     &,				(toev/toang)*Ez(iobs)
     	    write(iout,580)	twosigma
     &,				twosigma
     &,				twosigma
            write(iout,590)
     &				EFGxx(iobs)
     &,				EFGxy(iobs)
     &,				EFGxz(iobs)
     &,				EFGyy(iobs)
     &,				EFGyz(iobs)
     &,				EFGzz(iobs)
  560	    format(2x,'induced phi(r/A)/Volt = '
     &             ,	/'(+) ',e11.4,' (VOLT) +/- 'e11.4)
  565	    format(2x,
     $             'induced phi(r/A)/Volt due to surface charges = '
     &             ,	/'(+) ',e11.4,' (VOLT) +/- 'e11.4)
  570	    format('(-) ',e11.4,' * x',4x
     &             ,	'(-) ',e11.4,' * y',6x,'(-) ',e11.4,' * z')
  580	    format('+/- ',e11.4,' * x',4x
     &             ,	'+/- ',e11.4,' * y',6x,'+/- ',e11.4,' * z')
  590	    format('(-) ',e11.4,' * x^2',2x,'(-) ',e11.4,' * x * y',2x
     &             ,'(-) ',e11.4,' * x * z'
     &             ,	/23x,'(-) ',e11.4,' * y^2',4x,'(-) ',e11.4,
     &             ' *y * z',/48x,'(-) ',e11.4,' * z^2')
  160    continue
      endif
c
c
      if(prnt) then
         write(iout,1000) tokcal*deltag2
         write(iout,1010) tokcal*pomf2
         write(iout,1020) tokcal*aben
         write(iout,1030) tokcal*self
      endif
      if(debug) then
         do 124 i=1,nobs
            write(iout,'(4f10.5,2x)') q(i),qd(i,1),qd(i,2),qd(i,3)
            write(iout,'(6f10.5,2x)')  
     &                   qq(i,1),qq(i,2),qq(i,3),qq(i,4),qq(i,5),qq(i,6)
  124    continue
      endif
c
c
      return
      end
