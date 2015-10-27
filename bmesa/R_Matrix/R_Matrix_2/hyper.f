*deck hyper.f
c***begin prologue     hyper
c***date written       951229   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           
c***author             schneider, barry (nsf)
c***source             
c***purpose            r matrix scattering code for 2 dimensions.
c***                   
c***references         
c
c***routines called    
c***end prologue       hyper
      subroutine hyper(pham,energy,rbox,sym,n,nen,smatrix,dim,prn)
      implicit integer (a-z)
      integer*8 pham, ph, pgam, phamx0, phamy0
      integer*8 prmat, pkmat, psol, pscr
      real*8 rmat, ham, gam, hamx0, hamy0
      real*8 kmat, sol, smat, scr, rbox, energy
      character*80 title
      character*(*) sym
      character*16 fptoc
      logical prn, smatrix
      dimension pham(4,2), n(4,2), energy(nen), prn(*), ngot(6)
      common/io/inp, iout
      pointer (ph,ham(1))
      pointer (pgam,gam(1))
      pointer (phamx0,hamx0(1))
      pointer (phamy0,hamy0(1))
      pointer (prmat,rmat(1))
      pointer (pkmat,kmat(1))
      pointer (psol,sol(1))
      pointer (psol,isol(1))
      pointer (pscr,scr(1))
      pointer (pscr,iscr(1)) 
      pointer (psmat,smat(1))
c
c     calculate the r-matrix amplitudes
c
      call fgam(pham,pgam,sym,ngot(1),n,dim,prn(1))
      call lnkerr('quit')
      phamx0=pham(1,2)
      phamy0=pham(2,2)
      ph=pham(dim+1,1)
      hmtx0=1
      vx0=hmtx0+n(1,2)*n(1,2)
      eigvcx0=vx0+n(1,2)
      gmx0=eigvcx0+n(1,2)*n(1,2)      
      eigvlx0=gmx0+n(1,2)
      srfx0=eigvlx0+n(1,2)
      hmty0=1
      vy0=hmty0+n(2,2)*n(2,2)
      eigvcy0=vy0+n(2,2)
      gmy0=eigvcy0+n(2,2)*n(2,2)      
      eigvly0=gmy0+n(2,2)
      srfy0=eigvly0+n(2,2)
      eigvc=1
      eigvl=eigvc+n(dim+1,1)*n(dim+1,1)
c
c      the total energy must be larger than the lowest energy of the initial
c      state
c
      if(sym.eq.'unsymmetric') then
c
c        we dimension rmat like an nc*nc matrix and then pass the
c        pieces.  same for the solutions.  thus when we are eventually
c        done we have one large r-matrix and can match simply.
c
         nc=n(1,2)+n(2,2)
         rmt=1
         rmtxx=1
         rmtxy=rmtxx+n(1,2)*nc
         rmtyx=rmtxx+n(1,2)
         rmtyy=rmtxy+n(1,2)
         echn=rmt+nc*nc 
         echnx=echn
         echny=echnx+n(1,2)
         need=wpadti(echn+nc)
         call getmem(need,prmat,ngot(2),'rmat',0)
         kmt=1
         need=wpadti(kmt+nc*nc)
         call getmem(need,pkmat,ngot(3),'kmat',0)
         sn=1
         snx=sn
         sny=snx+n(1,2)
         dsn=sn+nc
         dsnx=dsn
         dsny=dsnx+n(1,2) 
         cn=dsn+nc
         cnx=cn
         cny=cnx+n(1,2)
         dcn=cn+nc
         dcnx=dcn
         dcny=dcnx+n(1,2)
         typ=wpadti(dcn+nc)
         typx=typ
         typy=typx+n(1,2)  
         need=typ+nc
         call getmem(need,psol,ngot(4),'solution',0)
         fac=1
         if(smatrix) then
            fac=2
            smt=1
            need=wpadti(smt+fac*nc*nc)
            call getmem(need,psmat,ngot(5),'smat',0)
         endif
         rc=1
         work=rc+fac*nc*nc
         lwork=5*fac*nc
         ipvt=wpadti(work+lwork)
         need=ipvt+nc
         call getmem(need,pscr,ngot(6),'scr',0)      
         gma=1
         gammax=gma
         gammay=gammax+n(1,2)
         do 20 ene=1,nen
            no=0
c
c           compute the r-matrices
c
            call rzero(rmat(rmt),nc*nc)
            call conrmat(rmat(rmtxx),rmat(rmtxx),ham(eigvl),
     1                   gam(gammax),gam(gammax),energy(ene),
     2                   n(1,2),n(1,2),nc,n(dim+1,1),.false.,.false.)
            call conrmat(rmat(rmtxy),rmat(rmtxy),ham(eigvl),
     1                   gam(gammax),gam(gammay),energy(ene),
     2                   n(1,2),n(2,2),nc,n(dim+1,1),.false.,.false.)
            call conrmat(rmat(rmtyx),rmat(rmtyx),ham(eigvl),
     1                   gam(gammay),gam(gammax),energy(ene),
     2                   n(2,2),n(1,2),nc,n(dim+1,1),.false.,.false.)
            call conrmat(rmat(rmtyy),rmat(rmtyy),ham(eigvl),
     1                   gam(gammay),gam(gammay),energy(ene),
     2                   n(2,2),n(2,2),nc,n(dim+1,1),.false.,.false.)
            if(prn(2)) then
               title='r-matrix gamma matrix'
               call prntrm(title,gam(gma),nc,n,nc,n,iout)
               title='R-matrix: energy = '//fptoc(energy(ene))
               call prntrm(title,rmat(rmt),nc,nc,nc,nc,iout)
            endif
c
c           compute the linearly independent solutions at the box
c
            call extrnl(sol(snx),sol(dsnx),sol(cnx),sol(dcnx),
     1                  hamx0(eigvlx0),rmat(echnx),isol(typx),rbox,
     2                  energy(ene),n(1,2),no,prn(3))
            call extrnl(sol(sny),sol(dsny),sol(cny),sol(dcny),
     1                  hamy0(eigvly0),rmat(echny),isol(typy),rbox,
     2                  energy(ene),n(2,2),no,prn(3))
            call kmtrx(rmat(rmt),sol(sn),sol(dsn),sol(cn),sol(dcn),
     1                 rmat(echn),isol(typ),energy(ene),kmat(kmt),
     2                 scr(rc),scr(work),iscr(ipvt),nc,no,lwork,prn(4))
            if(smatrix) then
               call smtrx(smat(smt),kmat(kmt),energy(ene),scr(rc),
     1                    scr(work),iscr(ipvt),nc,no,lwork,prn(5))
            endif
 20      continue   
      elseif(sym.eq.'symmetric') then
         nc=n(1,2)
         rmt=1
         echn=rmt+nc*nc
         need=wpadti(echn+nc)
         call getmem(need,prmat,ngot(2),'rmat',0)
         kmt=1
         need=wpadti(kmt+nc*nc)
         call getmem(need,pkmat,ngot(3),'kmat',0)
         sn=1
         dsn=sn+nc
         cn=dsn+nc
         dcn=cn+nc
         typ=wpadti(dcn+nc)
         need=wpadti(typ+n(1,2))
         call getmem(need,psol,ngot(4),'solution',0)
         fac=1
         if(smatrix) then
            fac=2
            smt=1
            need=wpadti(smt+fac*nc*nc)
            call getmem(need,psmat,ngot(5),'smat',0)
         endif
         rc=1
         work=rc+nc*nc*fac
         lwork=5*nc*fac
         ipvt=wpadti(work+lwork)
         need=ipvt+nc
         call getmem(need,pscr,ngot(6),'scr',0)      
         gamma=1
         no=0
         do 30 ene=1,nen
            call rzero(rmat(rmt),nc*nc)
            call conrmat(rmat(rmt),ham(eigvl),
     1                   gam(gamma),gam(gamma),
     2                   energy(ene),nc,nc,nc,n(dim+1,1),prn(2))
            call extrnl(sol(sn),sol(dsn),sol(cn),sol(dcn),
     1                  hamx0(eigvlx0),rmat(echn),isol(typ),rbox,
     2                  energy(ene),nc,no,prn(3))
            call kmtrx(rmat(rmt),sol(sn),sol(dsn),sol(cn),sol(dcn),
     1                 rmat(echn),isol(typ),energy(ene),
     2                 kmat(kmt),scr(rc),scr(work),iscr(ipvt),
     3                 nc,no,lwork,prn(4))
            if(smatrix) then
               call smtrx(smat(smt),kmat(kmt),energy(ene),scr(rc),
     1                    scr(work),iscr(ipvt),nc,no,lwork,prn(5))
            endif
 30      continue   
      else
         call lnkerr('error in symmetry type')
      endif
      call getmem(-ngot(1),pgam,idum,'gamma',idum) 
      call getmem(-ngot(2),prmat,idum,'rmat',idum)
      call getmem(-ngot(3),pkmat,idum,'kmat',idum)
      call getmem(-ngot(4),psol,idum,'solution',idum)
      call getmem(-ngot(6),pscr,idum,'scratch',idum)
      if(smatrix) then
         call getmem(-ngot(5),psmat,idum,'smat',idum)
      endif
      return      
      end       






