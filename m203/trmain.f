*deck @(#)trmain.f	1.2  7/30/91
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
      common/io/inp,iout
c
      write(iout,3000)
 3000 format(/,' m203: ',//,
     $' transform force constants from cartesian to ',
     $ 'internal coordinates')
c
      nat3=3*natoms
      step=0.005d+00
c
      call rzmat('rwf',nz,nvar,ianz,iz,bl,alpha,beta,
     $            klbl,lalpha,lbeta)
c
c
      call iosys('read real zvalues from rwf',nvar,xx(1,1),0,' ')
c
c     get the cartesian coordinates of the reference point
c
c
      call subvar(bl,alpha,beta,klbl,
     $            lalpha,lbeta,xx(1,1),nz,nvar)
c
c
      call ztoc(nz,ianz,iz,bl,alpha,beta,.false.,numtet,
     $          natoms,ian,atmchg,cref,scr6,scr1,scr2,scr3,
     $          scr4,scr5)
c
      write(iout,1050)
 1050 format(/' reference cartesian geometry'/)
      do 44 i=1,natoms
         k=3*(i-1)+1
         kk=k+2
         write(iout,2050) (cref(j),j=k,kk)
 2050    format(10x,3f18.5)
   44 continue
c
c           form the cartesian-to-internal transformation matrix
c
c
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
c
      call iosys('read real "cartesian second derivatives" from rwf',
     $            nvvc,ftric,0,' ')
c
c     write(iout,1992)
 1992 format(/'  b matrix ')
c      call matout(b,3*natoms,nvar,3*natoms,nvar,iout)
c
c     gradient independent transformation
c
c
      call trtosq(fsqc,ftric,3*natoms,nvvc)
c
      write(iout,1994)
 1994 format(/'  cartesian force constant matrix ')
      call matout(fsqc,3*natoms,3*natoms,3*natoms,3*natoms,iout)
c
      call ebc(scr11,fsqc,b,3*natoms,3*natoms,nvar)
      call ebtc(fsqi,b,scr11,nvar,3*natoms,nvar)
      write(iout,1995)
 1995 format(/' internal coordinate force constant matrix',/
     $,       '  gradient independent part ')
      call matout(fsqi,nvar,nvar,nvar,nvar,iout)
c
c
c     -- now calculate gradient dependent part of the transformation --
c     -- form a derivative b-matrix for each internal coordinate      --
c
      call iosys('read real "cartesian first derivatives" from rwf',
     $            nat3,gc,0,' ')
c
c      write(iout,1254) (gc(i),i=1,nat3)
 1254 format(' cartesian gradient ',/,9f8.4)
c
      do 5000 ivar=1,nvar
c
         xx(ivar,1)=xx(ivar,1) + step
c
         call cztran(nz,nvar,ianz,iz,bl,
     $            alpha,beta,tetrat,numtet,natoms,
     $            klbl,lalpha,lbeta,ian,atmchg,c,scr6,
     $            scr1,scr2,scr3,scr4,scr5,bplus,xx,cref,
     $            scr9,scr10)
c
         xx(ivar,1)=xx(ivar,1)- 2*step
c
         call cztran(nz,nvar,ianz,iz,bl,
     $            alpha,beta,tetrat,numtet,natoms,
     $            klbl,lalpha,lbeta,ian,atmchg,c,scr6,
     $            scr1,scr2,scr3,scr4,scr5,bminus,xx,cref,
     $            scr9,scr10)
c
         xx(ivar,1)=xx(ivar,1) + step
c
         do 5050 jj=1,nvar
            do 5040 ii=1,nat3
               bder(ii,jj)=(bplus(ii,jj)-bminus(ii,jj))/(2*step)
 5040       continue
 5050    continue
c
c
cc      write(iout,1999) ivar
 1999 format(/'  bplus,bminus,bder: ivar =  ',i3)
c      call matout(bplus,nat3,nvar,nat3,nvar,iout)
c      call matout(bminus,nat3,nvar,nat3,nvar,iout)
c      call matout(bder,nat3,nvar,nat3,nvar,iout)
c
         call ebtc(tt,bder,gc,nvar,nat3,1)
c      write(iout,1253) (tt(i),i=1,nvar)
 1253 format(/' tt ',/,3f20.7)
         do 5060 ll=1,nvar
            fsqi(ll,ivar)=fsqi(ll,ivar)+tt(ll)
 5060    continue
c
 5000 continue
      write(iout,1996)
 1996 format(/' internal coordinate force constant matrix',/)
      call matout(fsqi,nvar,nvar,nvar,nvar,iout)
c
      call sqtotr(frcnst,fsqi,nvar,nvv)
      call iosys('write real force_constants to rwf',nvv,frcnst,0,' ')
c
      return
      end
