*deck @(#)cztran.f	5.1  11/6/94
      subroutine cztran(nz,nvar,ianz,iz,bl,
     $             alpha,beta,numtet,natoms,
     $             klbl,lalpha,lbeta,ian,atmchg,c,scr6,
     $             scr1,scr2,scr3,scr4,scr5,b,xx,cref,
     $             cplus,cminus)
c***begin prologue     cztran.f
c***date written       yymmdd  
c***revision date      11/6/94      
c
c***keywords           
c***author             page, michael(nrl)
c***source             @(#)cztran.f	5.1   11/6/94
c***purpose            transform cartesian force constants to internals.
c***description
c     
c    
c
c***references
c
c***routines called
c
c***end prologue       cztran.f
      implicit none
c     --- input variables -----
      integer nz,nvar,numtet,natoms
c     --- input arrays (unmodified) ---
      integer lalpha(nz),lbeta(nz),klbl(nz)
      integer ianz(nz),iz(4,nz),ian(nz)
      real*8 bl(nz),alpha(nz),beta(nz),c(natoms*3),atmchg(nz)
c     --- input arrays (scratch) ---
      real*8 scr1(nz),scr2(nz),scr3(nz),scr4(nz),scr5(nz),scr6(nz*3)
      real*8 cplus(3*natoms),cminus(3*natoms)
c     --- output arrays ---
      real*8 b(3*natoms,nvar),xx(nvar,2*nvar+1),cref(3*natoms)
c     --- output variables ---
c     --- scratch arrays ---
c     --- local variables ---
      integer inp,iout
      integer ii,jj,iplus,iminus,k
      real*8 step,tenm3,two
      parameter (tenm3=0.001d+00,two=2.0d+00)
c
      common/io/inp,iout
c
c     --- loop over internal coordinates and generate cartesian
c         displacement vectors associated with small internal
c         coordinate displacements
c
c     --- first generate internal coordinate displacement vectors
c         for differencing
      step=tenm3
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
c     --- double point
      do 40 jj=2,nvar+1
         iplus=2*(jj-1)
         iminus=iplus+1
         call subvar(bl,alpha,beta,klbl,
     $               lalpha,lbeta,xx(1,iplus),nz,nvar)
         call ztoc(nz,ianz,iz,bl,alpha,beta,.false.,numtet,
     $             natoms,ian,atmchg,cplus,scr6,scr1,scr2,scr3,scr4,
     $             scr5)
         call subvar(bl,alpha,beta,klbl,
     $               lalpha,lbeta,xx(1,iminus),nz,nvar)
         call ztoc(nz,ianz,iz,bl,alpha,beta,.false.,numtet,
     $             natoms,ian,atmchg,cminus,scr6,scr1,scr2,scr3,scr4,
     $             scr5)
c
         do 60 k=1,3*natoms
            b(k,jj-1)=(cplus(k)-cminus(k))/(two*step)
   60    continue
   40 continue
c
c     --- b contains the derivatives of cartesian coordinates 
c         with respect to internals.
c
      return
      end
