*deck @(#)trmain.f	5.1  11/6/94
      subroutine trmain(nvar,nvv,nz,natoms,nvvc,toang,ops,x,f,c,
     $                fsqc,ftric,xx,frcnst,vname,lbl,ian,ianz,
     $                iz,bl,alpha,beta,klbl,lalpha,lbeta,scr1,scr2,
     $                scr3,scr4,scr5,scr6,bplus,bminus,bder,tt,gc,
     $                b,cref,scr9,scr10,scr11,fsqi,atmchg)
c***begin prologue     trmain.f
c***date written       yymmdd  
c***revision date      11/6/94      
c
c***keywords           
c***author             page, michael(nrl)
c***source             @(#)trmain.f	5.1   11/6/94
c***purpose            transform force constants from cartesian to internal.
c***description
c     
c    
c
c***references
c
c***routines called
c
c***end prologue       trmain.f
      implicit none
c     --- input variables -----
      integer nvar,nz,nvv,nvvc,natoms
      character*(*) ops
c     --- input arrays (unmodified) ---
      character*(*) vname(nvar)
      integer lbl(nz),lalpha(nz),lbeta(nz),klbl(nz)
      integer ianz(nz),iz(4,nz),ian(nz)
c     --- input arrays (scratch) ---
c     --- output arrays ---
      real*8 x(nvar),f(nvar),frcnst(nvv)
      real*8 ftric(nvvc),fsqc(3*natoms,3*natoms)
      real*8 xx(nvar,2*nvar+1),fsqi(nvar,nvar)
      real*8 bl(nz),alpha(nz),beta(nz),atmchg(nz)
      real*8 c(3*natoms),tt(nvar),gc(3*natoms)
      real*8 scr1(nz),scr2(nz),scr3(nz),scr4(nz),scr5(nz),scr6(3*nz)
      real*8 b(3*natoms,nvar),cref(3*natoms),scr9(3*natoms)
      real*8 scr10(3*natoms),scr11(nvar*3*natoms)
      real*8 bplus(3*natoms,nvar),bminus(3*natoms,nvar)
      real*8 bder(3*natoms,nvar)
c     --- output variables ---
c     --- scratch arrays ---
c     --- local variables ---
      integer inp,iout
      integer nat3,nvc,i,j,k,kk,ii,jj,ll,ivar
      integer numtet
      real*8 toang,step,hftenm3
      logical debug
      parameter (debug=.false.,hftenm3=0.005d+00)
c
      common/io/inp,iout
c
 1000 format(1x,'m732:',
     $          'transform force constants from cartesian to ',
     $          'internal coordinates')
 1050 format(5x,'reference cartesian geometry')
 1100 format(10x,3f18.5)
 1150 format(/'  b matrix ')
 1200 format(5x,'cartesian force constant matrix ')
 1250 format(5x,'internal coordinate force constant matrix:',
     $          'gradient independent part ')
 1300 format(5x,'cartesian gradient',/,9f8.4)
 1350 format(/'  bplus,bminus,bder: ivar =  ',i3)
 1400 format(/' tt ',/,3f20.7)
 1450 format(5x,'internal coordinate force constant matrix')
c
c     --- announce our presence
      write(iout,1000)
c
c     --- let's go, retrieve z-matrix variables
      nat3=3*natoms
      step=hftenm3
      call rzmat('rwf',nz,nvar,ianz,iz,bl,alpha,beta,
     $            klbl,lalpha,lbeta)
      call iosys('read real zvalues from rwf',nvar,xx(1,1),0,' ')
c
c     --- get the cartesian coordinates of the reference point
      call subvar(bl,alpha,beta,klbl,
     $            lalpha,lbeta,xx(1,1),nz,nvar)
      call ztoc(nz,ianz,iz,bl,alpha,beta,.false.,numtet,
     $          natoms,ian,atmchg,cref,scr6,scr1,scr2,scr3,
     $          scr4,scr5)
      if(debug) then
         write(iout,1050)
         do 10 i=1,natoms
            k=3*(i-1)+1
            kk=k+2
            write(iout,1100) (cref(j),j=k,kk)
   10    continue
      endif
c
c     --- form the cartesian-to-internal transformation matrix
      call cztran(nz,nvar,ianz,iz,bl,
     $            alpha,beta,numtet,natoms,
     $            klbl,lalpha,lbeta,ian,atmchg,c,scr6,
     $            scr1,scr2,scr3,scr4,scr5,b,xx,cref,
     $            scr9,scr10)
c
      nvc=nvar*3*natoms
      call iosys('write real "ctoz transformation matrix" to rwf',
     $            nvc,b,0,' ')
c
c    --- read the cartesian force constant matrix from rwf ---
      call iosys('read real "cartesian second derivatives" from rwf',
     $            nvvc,ftric,0,' ')
      if(debug) then
         write(iout,1150)
         call matout(b,3*natoms,nvar,3*natoms,nvar,iout)
      endif
c
c     --- gradient independent transformation
      call trtosq(fsqc,ftric,3*natoms,nvvc)
      if(debug) then
         write(iout,1200)
         call matout(fsqc,3*natoms,3*natoms,3*natoms,3*natoms,iout)
      endif
      call ebc(scr11,fsqc,b,3*natoms,3*natoms,nvar)
      call ebtc(fsqi,b,scr11,nvar,3*natoms,nvar)
      if(debug) then
         write(iout,1250)
         call matout(fsqi,nvar,nvar,nvar,nvar,iout)
      endif
c
c     --- now calculate gradient dependent part of the transformation --
c         form a derivative b-matrix for each internal coordinate      --
      call iosys('read real "cartesian first derivatives" from rwf',
     $            nat3,gc,0,' ')
      if(debug) then
         write(iout,1300) (gc(i),i=1,nat3)
      endif
      do 100 ivar=1,nvar
         xx(ivar,1)=xx(ivar,1) + step
         call cztran(nz,nvar,ianz,iz,bl,
     $            alpha,beta,numtet,natoms,
     $            klbl,lalpha,lbeta,ian,atmchg,c,scr6,
     $            scr1,scr2,scr3,scr4,scr5,bplus,xx,cref,
     $            scr9,scr10)
         xx(ivar,1)=xx(ivar,1)- 2*step
         call cztran(nz,nvar,ianz,iz,bl,
     $            alpha,beta,numtet,natoms,
     $            klbl,lalpha,lbeta,ian,atmchg,c,scr6,
     $            scr1,scr2,scr3,scr4,scr5,bminus,xx,cref,
     $            scr9,scr10)
         xx(ivar,1)=xx(ivar,1) + step

         do 50 jj=1,nvar
            do 40 ii=1,nat3
               bder(ii,jj)=(bplus(ii,jj)-bminus(ii,jj))/(2*step)
 40         continue
 50      continue
c
         if(debug) then
            write(iout,1350) ivar
            call matout(bplus,nat3,nvar,nat3,nvar,iout)
            call matout(bminus,nat3,nvar,nat3,nvar,iout)
            call matout(bder,nat3,nvar,nat3,nvar,iout)
         endif
         call ebtc(tt,bder,gc,nvar,nat3,1)
c
         if(debug) then
            write(iout,1400) (tt(i),i=1,nvar)
         endif
c
         do 60 ll=1,nvar
            fsqi(ll,ivar)=fsqi(ll,ivar)+tt(ll)
   60    continue
c
  100 continue
c
c
      if(debug) then
         write(iout,1450)
         call matout(fsqi,nvar,nvar,nvar,nvar,iout)
      endif
      call sqtotr(frcnst,fsqi,nvar,nvv)
      call iosys('write real force_constants to rwf',nvv,frcnst,0,' ')
c
c
      return
      end
