*deck %W%  %G%
      subroutine trmain(nvar,nvv,nz,natoms,nvvc,toang,ops,x,f,c,
     $                fsqc,ftric,xx,frcnst,vname,lbl,ian,ianz,
     $                iz,bl,alpha,beta,klbl,lalpha,lbeta,scr1,scr2,
     $                scr3,scr4,scr5,scr6,bplus,bminus,bder,tt,gc,
     $                b,cref,scr9,scr10,scr11,fsqi,atmchg)
c
      implicit integer(a-z)
      character*(*) ops,vname(nvar)
      real*8 x(nvar),f(nvar),frcnst(nvv)
      real*8 ftric(nvvc),fsqc(3*natoms,3*natoms)
      real*8 xx(nvar,2*nvar+1),fsqi(nvar,nvar)
      real*8 toang,step
      real*8 bl(nz),alpha(nz),beta(nz),atmchg(nz)
      real*8 c(3*natoms),tt(nvar),gc(3*natoms)
      real*8 scr1(nz),scr2(nz),scr3(nz),scr4(nz),scr5(nz),scr6(3*nz)
      real*8 b(3*natoms,nvar),cref(3*natoms),scr9(3*natoms)
      real*8 scr10(3*natoms),scr11(nvar*3*natoms)
      real*8 bplus(3*natoms,nvar),bminus(3*natoms,nvar)
      real*8 bder(3*natoms,nvar)
      integer lbl(nz),lalpha(nz),lbeta(nz),klbl(nz)
      integer ianz(nz),iz(4,nz),ian(nz)
      logical debug
      parameter (debug=.false.)
      common/io/inp,iout
c
      write(iout,3000)
 3000 format(1x,'m203:',
     $          'transform force constants from cartesian to ',
     $          'internal coordinates')
c
c     --- let's go, retrieve z-matrix variables
      nat3=3*natoms
      step=0.005d+00
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
 1050    format(5x,'reference cartesian geometry')
         do 44 i=1,natoms
            k=3*(i-1)+1
            kk=k+2
            write(iout,2050) (cref(j),j=k,kk)
 2050       format(10x,3f18.5)
   44    continue
      endif
c
c     --- form the cartesian-to-internal transformation matrix
      call cztran(nz,nvar,ianz,iz,bl,
     $            alpha,beta,tetrat,numtet,natoms,
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
         write(iout,1992)
 1992    format(/'  b matrix ')
         call matout(b,3*natoms,nvar,3*natoms,nvar,iout)
      endif
c
c     --- gradient independent transformation
      call trtosq(fsqc,ftric,3*natoms,nvvc)
      if(debug) then
         write(iout,1994)
 1994    format(5x,'cartesian force constant matrix ')
         call matout(fsqc,3*natoms,3*natoms,3*natoms,3*natoms,iout)
      endif
      call ebc(scr11,fsqc,b,3*natoms,3*natoms,nvar)
      call ebtc(fsqi,b,scr11,nvar,3*natoms,nvar)
      if(debug) then
         write(iout,1995)
 1995    format(5x,'internal coordinate force constant matrix:',
     $             'gradient independent part ')
         call matout(fsqi,nvar,nvar,nvar,nvar,iout)
      endif
c
c     --- now calculate gradient dependent part of the transformation --
c         form a derivative b-matrix for each internal coordinate      --
      call iosys('read real "cartesian first derivatives" from rwf',
     $            nat3,gc,0,' ')
      if(debug) then
         write(iout,1254) (gc(i),i=1,nat3)
 1254    format(5x,'cartesian gradient',/,9f8.4)
      endif
      do 5000 ivar=1,nvar
         xx(ivar,1)=xx(ivar,1) + step
         call cztran(nz,nvar,ianz,iz,bl,
     $            alpha,beta,tetrat,numtet,natoms,
     $            klbl,lalpha,lbeta,ian,atmchg,c,scr6,
     $            scr1,scr2,scr3,scr4,scr5,bplus,xx,cref,
     $            scr9,scr10)
         xx(ivar,1)=xx(ivar,1)- 2*step
         call cztran(nz,nvar,ianz,iz,bl,
     $            alpha,beta,tetrat,numtet,natoms,
     $            klbl,lalpha,lbeta,ian,atmchg,c,scr6,
     $            scr1,scr2,scr3,scr4,scr5,bminus,xx,cref,
     $            scr9,scr10)
         xx(ivar,1)=xx(ivar,1) + step

         do 5050 jj=1,nvar
            do 5040 ii=1,nat3
               bder(ii,jj)=(bplus(ii,jj)-bminus(ii,jj))/(2*step)
 5040       continue
 5050    continue
c
         if(debug) then
            write(iout,1999) ivar
 1999       format(/'  bplus,bminus,bder: ivar =  ',i3)
            call matout(bplus,nat3,nvar,nat3,nvar,iout)
            call matout(bminus,nat3,nvar,nat3,nvar,iout)
            call matout(bder,nat3,nvar,nat3,nvar,iout)
         endif
         call ebtc(tt,bder,gc,nvar,nat3,1)
c
         if(debug) then
            write(iout,1253) (tt(i),i=1,nvar)
 1253       format(/' tt ',/,3f20.7)
         endif
c
         do 5060 ll=1,nvar
            fsqi(ll,ivar)=fsqi(ll,ivar)+tt(ll)
 5060    continue
c
 5000 continue
c
c
      if(debug) then
         write(iout,1996)
 1996    format(5x,'internal coordinate force constant matrix')
         call matout(fsqi,nvar,nvar,nvar,nvar,iout)
      endif
      call sqtotr(frcnst,fsqi,nvar,nvv)
      call iosys('write real force_constants to rwf',nvv,frcnst,0,' ')
c
c
      return
      end
