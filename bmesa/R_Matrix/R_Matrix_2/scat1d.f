*deck scat1d.f
c***begin prologue     scat1d
c***date written       951229   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           
c***author             schneider, barry (nsf)
c***source             
c***purpose            one dimensional scattering.
c***                   
c***references         
c
c***routines called    
c***end prologue       scat1
      subroutine scat1d(pham,energy,rbox,n,nen,smatrix,prn)
      implicit integer (a-z)
      integer*8 pham, prmat, pkmat, psol, pscr, psmat 
      real*8 rmat, ham, kmat, sol, scr, energy, rbox, pse, pi
      real*8 sinpse, xsect
      complex*16 smat
      logical prn, smatrix
      dimension energy(nen), prn(3), ngot(5)
      common/io/inp, iout
      data pi  / 3.141592653589793238462643d+00 /
      pointer (pham,ham(1))
      pointer (prmat,rmat(1))
      pointer (pkmat,kmat(1))
      pointer (psol,sol(1))
      pointer (psol,isol(1))
      pointer (pscr,scr(1))
      pointer (pscr,iscr(1)) 
      pointer (psmat,smat(1))
c
      nc=1
      no=nc
      rmt=1
      echn=rmt+nc*nc
      need=wpadti(echn+nc)
      call getmem(need,prmat,ngot(1),'rmat',0)         
      kmt=1
      need=wpadti(kmt+nc*nc)
      call getmem(need,pkmat,ngot(2),'kmat',0)         
      sn=1
      dsn=sn+nc
      cn=dsn+nc
      dcn=cn+nc
      typ=wpadti(dcn+nc)
      need=typ+nc
      call getmem(need,psol,ngot(3),'solution',0)
      fac=1
      if(smatrix) then
         fac=2
         smt=1
         need=wpadti(smt+fac*no*no)
         call getmem(need,psmat,ngot(4),'smat',0)
      endif
      rc=1
      work=rc+fac*nc*nc
      lwork=5*fac*nc
      ipvt=wpadti(work+lwork)
      need=ipvt+nc
      call getmem(need,pscr,ngot(5),'scratch',0)
      lwork=lwork/fac
c
c     set the pointers to the r-matrix amplitudes.  in one dimension
c     these are already available from the call to lobatto.
c
      eigval=1
      gamma=eigval+n
      do 10 ene=1,nen
         call rzero(rmat(rmt),nc*nc)
         call conrmat(rmat(rmt),rmat(rmt),ham(eigval),
     1                ham(gamma),ham(gamma),energy(ene),
     2                nc,nc,nc,n,prn(1),.false.)
         no=0
         call extrnl(sol(sn),sol(dsn),sol(cn),sol(dcn),0.d0,rmat(echn),
     1               isol(typ),rbox,energy(ene),nc,no,prn(3))
         call kmtrx(rmat(rmt),sol(sn),sol(dsn),sol(cn),sol(dcn),
     1              rmat(echn),isol(typ),energy(ene),
     2              kmat(kmt),scr(rc),scr(work),iscr(ipvt),nc,
     3              no,lwork,prn(4))
         pse=atan(kmat(kmt))
         sinpse=sin(pse)
         xsect=2.d0*pi*sinpse*sinpse/(energy(ene))
         write(iout,1) energy(ene), pse, xsect
         if(smatrix) then
            call smtrx1(smat(smt),rmat(rmt),energy(ene),rbox)
         endif
 10   continue   
      call getmem(-ngot(1),prmat,idum,'rmat',idum)
      call getmem(-ngot(2),pkmat,idum,'kmat',idum)
      call getmem(-ngot(3),psol,idum,'solution',idum)
      call getmem(-ngot(5),pscr,idum,'scratch',idum)
      if(smatrix) then
         call getmem(-ngot(4),psmat,idum,'smat',idum)
      endif
      return      
 1    format(/,1x,'    energy        =     ',e15.8,
     1            '    phase shift   =     ',e15.8,/,1x,
     2            '    cross section =     ',e15.8)
      end       


















