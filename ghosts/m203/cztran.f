*deck %W%  %G%
c
      subroutine cztran(nz,nvar,ianz,iz,bl,
     $             alpha,beta,tetrat,numtet,natoms,
     $             klbl,lalpha,lbeta,ian,atmchg,c,scr6,
     $             scr1,scr2,scr3,scr4,scr5,b,xx,cref,
     $             cplus,cminus)
      implicit integer(a-z)
      real*8 bl(nz),alpha(nz),beta(nz),c(natoms*3),atmchg(nz)
      real*8 scr1(nz),scr2(nz),scr3(nz),scr4(nz),scr5(nz),scr6(nz*3)
      real*8 b(3*natoms,nvar),xx(nvar,2*nvar+1),cref(3*natoms)
      real*8 cplus(3*natoms),cminus(3*natoms)
      real*8 step
      dimension lalpha(nz),lbeta(nz),klbl(nz)
      dimension ianz(nz),iz(4,nz),ian(nz)
      common/io/inp,iout
c
c
c
c     loop over internal coordinates and generate cartesian
c     displacement vectors associated with small internal
c     coordinate displacements
c
c     first generate internal coordinate displacement vectors
c     for differencing
      step=.001d+00
c
      do 35 jj=2,nvar+1
         iplus=2*(jj-1)
         iminus=iplus+1
           do 36 ii=1,nvar
             xx(ii,iplus)=xx(ii,1)
             xx(ii,iminus)=xx(ii,1)
   36      continue
         xx(jj-1,iplus)=xx(jj-1,1)+step
         xx(jj-1,iminus)=xx(jj-1,1)-step
   35 continue
c
c
c     double point
c
      do 40 jj=2,nvar+1
c
c      write(iout,1200) (xx(i,jj),i=1,nvar)
            iplus=2*(jj-1)
            iminus=iplus+1
c
            call subvar(bl,alpha,beta,klbl,
     $                  lalpha,lbeta,xx(1,iplus),nz,nvar)
c
            call ztoc(nz,ianz,iz,bl,alpha,beta,.false.,numtet,
     $                natoms,ian,atmchg,cplus,scr6,scr1,scr2,scr3,scr4,
     $                scr5)
c
            call subvar(bl,alpha,beta,klbl,
     $                  lalpha,lbeta,xx(1,iminus),nz,nvar)
c
            call ztoc(nz,ianz,iz,bl,alpha,beta,.false.,numtet,
     $                natoms,ian,atmchg,cminus,scr6,scr1,scr2,scr3,scr4,
     $                scr5)
c
      do 60 k=1,3*natoms
         b(k,jj-1)=(cplus(k)-cminus(k))/(2.0*step)
   60 continue
c      write(iout,1775)
 1775 format(' cplus and cminus , ')
cc      write(iout,2050) (cplus(ll),ll=1,3*natoms)
cc      write(iout,2050) (cminus(ll),ll=1,3*natoms)
c
   40 continue
c
c      write(iout,1992)
 1992 format(/'  b matrix ')
c      call matout(b,3*natoms,nvar,3*natoms,nvar,iout)
c
c     -- b contains the derivatives of cartesian coordinates --
c     -- with respect to internals.
c
      return
      end
