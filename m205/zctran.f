*deck @(#)zctran.f	5.1  11/6/94
      subroutine zctran(nz,nvar,ianz,iz,bl,
     $             alpha,beta,tetrat,numtet,natoms,
     $             klbl,lalpha,lbeta,ian,atmchg,c,scr6,
     $             scr1,scr2,scr3,scr4,scr5,b,xx,cref,
     $             cplus,cminus,order,binv)
c***begin prologue     zctran.f
c***date written       yymmdd  
c***revision date      11/6/94      
c
c***keywords           
c***author             page, michael(nrl)
c***source             @(#)zctran.f	5.1   11/6/94
c***purpose            
c***description
c     
c    
c
c***references
c
c***routines called
c
c***end prologue       zctran.f
      implicit none
c     --- input variables -----
      integer nz,nvar,numtet,natoms
c     --- input arrays (unmodified) ---
      real*8 atmchg(nz),c(natoms*3)
c     --- input arrays (scratch) ---
c     --- output arrays ---
      real*8 b(3*natoms,nvar)
c     --- output variables ---
c     --- scratch arrays ---
      integer lalpha(nz),lbeta(nz),klbl(nz)
      integer ianz(nz),iz(4,nz),ian(nz),order(3*natoms)
      real*8 bl(nz),alpha(nz),beta(nz)
      real*8 scr1(nz),scr2(nz),scr3(nz),scr4(nz),scr5(nz),scr6(nz*3)
      real*8 xx(nvar,2*nvar+1),cref(3*natoms)
      real*8 cplus(3*natoms),cminus(3*natoms),binv(nvar,nvar)
c     --- local variables ---
      integer d2ecycl
      integer i,j,k,kk,jj,iplus,iminus
      integer inp,iout,irt
      logical prnt,chkpt,singpt,cartfx,tetrat
      logical debug
      real*8 det
      real*8 energy,rmax,rmin,rlim,stpsize
c
      parameter (debug=.false.)
c
      common/io/inp,iout
      common/d2einf/energy,rmax,rmin,rlim,d2ecycl,
     $               prnt,chkpt,singpt,stpsize,cartfx
c
 1000 format(5x,'cartesians '/,(3f14.4))
 1050 format(5x,'reference cartesian geometry'/)
 1200 format(5x,'internals: ',3f14.4)
 1555 format(5x,'c to z trans matrix - inverting'/)
 1556 format(5x,'after inversion')
 1559 format(5x,'b matrix before contraction '/)
 2050 format(10x,3f18.5)
c
c     --- retrieve z-matrix information
      call rzmat('rwf',nz,nvar,ianz,iz,bl,alpha,beta,
     $            klbl,lalpha,lbeta)
c
c     --- get the cartesian coordinates of the reference point
      call subvar(bl,alpha,beta,klbl,
     $            lalpha,lbeta,xx(1,1),nz,nvar)
c
      call ztoc(nz,ianz,iz,bl,alpha,beta,.false.,numtet,
     $          natoms,ian,atmchg,cref,scr6,scr1,scr2,scr3,
     $          scr4,scr5)
c
      if(debug) then
         write(iout,1050)
         do 44 i=1,natoms
            k=3*(i-1)+1
            kk=k+2
            write(iout,2050) (cref(j),j=k,kk)
   44    continue
      endif
c
c     --- loop over internal coordinates and generate cartesian
c     displacement vectors associated with small internal
c     coordinate displacements
      if(singpt) then
         do 10 jj=2,nvar+1
            if(debug) then
               write(iout,1200) (xx(i,jj),i=1,nvar)
            endif
            call subvar(bl,alpha,beta,klbl,
     $                  lalpha,lbeta,xx(1,jj),nz,nvar)
            call ztoc(nz,ianz,iz,bl,alpha,beta,tetrat,numtet,
     $                natoms,ian,atmchg,c,scr6,scr1,scr2,scr3,scr4,
     $                scr5)
            do 30 k=1,3*natoms
               b(k,jj-1)=(c(k)-cref(k))/stpsize
   30       continue
            if(debug) then
               write(iout,1000) (c(i),i=1,3*natoms)
            endif
   10    continue
         if(debug) then
            write(iout,1559)
            call matout(b,natoms*3,nvar,natoms*3,nvar,iout)
         endif
      else
c        --- double point
         do 40 jj=2,nvar+1
            if(debug) then
               write(iout,1200) (xx(i,jj),i=1,nvar)
            endif
            iplus=2*(jj-1)
            iminus=iplus+1
            call subvar(bl,alpha,beta,klbl,
     $                  lalpha,lbeta,xx(1,iplus),nz,nvar)
            call ztoc(nz,ianz,iz,bl,alpha,beta,.false.,numtet,
     $                natoms,ian,atmchg,cplus,scr6,scr1,scr2,scr3,scr4,
     $                scr5)
            call subvar(bl,alpha,beta,klbl,
     $                  lalpha,lbeta,xx(1,iminus),nz,nvar)
            call ztoc(nz,ianz,iz,bl,alpha,beta,.false.,numtet,
     $                natoms,ian,atmchg,cminus,scr6,scr1,scr2,scr3,scr4,
     $                scr5)
            do 60 k=1,3*natoms
               b(k,jj-1)=(cplus(k)-cminus(k))/(2.0*stpsize)
   60       continue
   40    continue
      endif
c
c     --- now b contains the derivatives of cartesian coordinates 
c         with respect to internals. we need to invert this, but 
c         first contract the rectangular matrix into a square.  
c         this removes the fixed cartesian coordinates
c
c         fill up the order array which allows us to refer to the fixed
c         coordinates (order(1)-order(6)) or to the varied coordinates
c         (order(7)-order(3n))
      do 27 i=1,3*natoms
         order(i)=i
   27 continue
      order(6)=8
      order(7)=6
      order(8)=7
c
      irt=6
      do 88 j=1,nvar
         jj=order(j+irt)
         do 87 i=1,nvar
           binv(j,i)=b(jj,i)
   87    continue
   88 continue
      if(debug) then
         write(iout,1555)
         call matout(binv,nvar,nvar,nvar,nvar,iout)
      endif
      call minvrt(binv,nvar,nvar,det,scr1,scr2)
      if(debug) then
         write(iout,1556)
         call matout(binv,nvar,nvar,nvar,nvar,iout)
      endif
c
c     --- now put it back
      do 86 j=1,nvar
         jj=order(j+irt)
         do 85 i=1,nvar
            b(jj,i)=binv(i,j)
   85    continue
   86 continue
c
c
      return
      end
