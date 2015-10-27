*deck %W% %G%
      subroutine ksham(s,t,v,h,smhalf,f,u,eigval,t1,t2,t3,triang,values,
     $     ian,coord,c,wts,nradial,nomega,natoms,nbf,nnp,d,dlast,
     $     t4,t5,error,focksv,mxiter,finish,stdiis,
     $     ntriang,bmatrx,nshell,
     $     shlmin,shlmax,alpha,beta,dodiis,pulay,page,calc,
     $     ncoul,nexch,ndmat,energs,diissv,prnt,ops,
     $     mxds,ptr,dsptr,level,fcoef,jmat,jlast,kmat,
     $     bflabl,usesym,symlabl,salc,numso,lambda,lirrep,nirrep,
     $     orbsym,occsym,temp,kept,left,mxgrd,ngrid,dmcut,
     $     dencut,ptprim,noprim,nbtype,
     $     ex,nx,ny,nz,lenxyz,
     $     nocart,mintyp,maxmom,mxcont,nprim,
     $     pstart,prtoao,dpr,jprim,npf,nnprim,
     $     dirj,poisj,ntotal,nonzer,
     $     cont,nocont,ptcont,ntypes,ncont,start,nobf,minmom,
     $     xyzgrid,charge,maxl,bigl,dograd,
     $     slater,becke,vwn,lyp,
     $     mxgbsiz,valuesi,names,grdtyp,
     $     vwts,rnuc,amu,pwtx,rr,radii,akl,ptrad,
     $     rmax,lmax,adjust,minesz,nnshl,qint,dijmax)
c***begin prologue     %M%
c***date written       930511   (yymmdd)   
c***revision date      %G%      
c
c***keywords           
c***author             martin, richard (lanl) 
c***source             %W% %G%
c***purpose            
c***description
c     
c    
c
c***references
c                      w.kohn and l.j.sham, phys.rev.a 140,1133(1965).
c
c***routines called
c
c***end prologue       %M%
      implicit none
c     --- input variables -----
      integer nradial,nomega,natoms,nbf,nnp,mxiter
      integer finish,stdiis
      integer ntriang,nshell,ncoul,nexch,ndmat,mxds
      integer nirrep,left,mxgrd
      integer ntypes,nbtype,lenxyz,nprim,ncont,mxcont
      integer npf,nnprim
      integer mxgbsiz,bigl
      integer rmax,lmax,nnshl
      logical adjust
      logical dodiis,pulay,prnt,page,dirj,dograd,poisj
      real*8 dencut,dmcut
c     --- input arrays (unmodified) ---
      integer ian(natoms)
      integer shlmin(nshell),shlmax(nshell),ngrid(natoms)
      integer numso(nirrep),lambda(nirrep),occsym(nirrep,nshell)
      integer ptprim(natoms,nbtype),noprim(natoms,nbtype)
      integer nocart(nbtype),mintyp(nbtype)
      integer ptcont(natoms,ntypes),start(natoms,ntypes)
      integer nobf(ntypes),minmom(ntypes)
      integer nocont(nbtype)
      integer maxmom(nbtype)
      integer nx(lenxyz),ny(lenxyz),nz(lenxyz)
      integer pstart(natoms,nbtype)
      character*(*) ops
      character*(*) bflabl(nbf)
      character*(*) symlabl(nbf)
      character*(*) lirrep(nirrep)
      character*(*) calc
      character*16 names(*),grdtyp(*)
      real*8 s(nnp),t(nnp),v(nnp),smhalf(nnp)
      real*8 coord(3,natoms)
      real*8 salc(nbf,nbf)
      real*8 alpha(nshell,nshell),beta(nshell,nshell),fcoef(nshell)
      real*8 ex(nprim),cont(ncont)
      real*8 prtoao(npf*nbf)
c     --- input arrays (scratch) ---
      integer ptr(mxiter),dsptr(mxds)
      integer orbsym(nbf),temp(nbf,2),kept(nshell)
      integer ntotal(mxiter),nonzer(mxiter)
      integer maxl(natoms)
      integer ptrad(rmax)
      integer minesz
      integer valuesi(*)
      real*8 h(nnp),u(nbf,nbf),eigval(nbf)
      real*8 t1(*),t2(nbf,nbf),t3(nbf,nbf),triang(nnp)
      real*8 values(left),c(nbf,nbf),bmatrx(*)
      real*8 t4(nbf,nbf),t5(nbf,nbf),error(nbf,nbf,mxds)
      real*8 focksv(nnp,mxds)
      real*8 f(nnp),d(nnp,ndmat),dlast(nnp,ndmat)
      real*8 jmat(nnp,ncoul),jlast(nnp,ncoul),kmat(nnp,nexch)
      real*8 energs(mxiter),diissv(mxiter)
      real*8 dpr(nnprim,ndmat)
      real*8 jprim(nnprim,ndmat)
      real*8 xyzgrid(mxgrd,3,natoms),wts(mxgrd,natoms)
      real*8 charge(natoms,ndmat)
      real*8 vwts(mxgrd),rnuc(natoms,natoms),amu(natoms,natoms)
      real*8 pwtx(natoms),rr(natoms),radii(natoms),akl(natoms,natoms)
      real*8 qint(nnshl),dijmax(nnshl)
c     --- output arrays ---
c     --- output variables ---
c     --- scratch arrays ---
c     --- local variables ---
      integer inp,iout,ierr,stderr
      integer i,j,ndiis,iter,mndiis,maxocc,shell,nprint
      integer tmpindx,itmp,frcolap
      integer intkey
      integer dtot,imaxl,itch,aitch
      integer wpadti,iadtwp
      integer pleft
      integer vorder,vlmax,vradial,vncrule
      integer nlm,rpts,vlm,y2
      character*16 chrkey
      logical logkey,printj,printk,first,usesym,debug,abort
      logical becke,slater,vwn,lyp,calce
      logical maxolap,prteexch
      logical rhotest
      logical timeit,old
      real*8 energy,enuc,eonel,etwoel,ecoul,exc,eexch
      real*8 zero,ten,two
      real*8 diiser,maxerr,level,elast,deltae
      real*8 damp(50)
      real*8 tmpmx,cutexp,fpkey,totchg,qtest
      real*8 timej(3),timek(3)
      real*8 dum1,dum2,dum3,dum4,dum5,dum6
      real*8 biggy
c
      common/io/inp,iout
c
      parameter (debug=.false.)
      parameter (old=.false.)
      parameter (zero=0.0d+00,two=2.0d+00,ten=1.0d+01)
      parameter (timeit=.false.)
c
      data damp/50*1./
      data first /.true./
      data maxerr/1.0d-06/
      save damp,first,maxerr
c
 1000 format(' damping will be invoked during the',
     $           ' kohn-sham iterations',/,
     $           ' damping factor per shell ',
     $           /,10(2x,f5.3))
 1010 format (' the kinetic-energy one-electron integrals')
 1020 format (' the potential-energy one-electron integrals')
 1025 format(5x,'direct formation of j-matrix'
     $     /,8x,'pre-exponential integral cutoff ',1pe9.1,
     $     /,8x,'integral estimator cutoff       ',1pe9.1)
 1026 format(5x,'directj needs more core')
 1027 format(5x,'have:',i8,' need:',i8, ' double precision words')
 1028 format(5x,'direct poisson solution for j-matrix',
     $      /8x,'maximum l:        ',i3,
     $      /8x,'nradial:          ',i3,
     $      /8x,'newton-cotes rule:',i3,
     $      /8x,'spline-order:     ',i3)
 1030 format (' the total one-electron integrals ')
 1040 format (' initial guess vector ')
 1050 format(5x,'nuclear repulsion energy:',7x,f12.6)
 1060 format(5x,'iter              energy    diis convergence')
 1065 format(5x,'iter              energy    diis convergence',
     $       3x,'integrals')
 1070 format (' density matrices')
 1080 format (' q-matrix in the m. o. basis')
 1090 format (' the fock-matrices in the m.o. basis ')
 1100 format (' the fock-matrices in the m.o. basis ')
 1110 format(5(2x,f10.4))
 1120 format (' eigenvectors and eigenvalues of f(tilde)')
 1130 format (' vector at the end of iteration',i5)
 1135 format(5x,i4,2(5x,f15.9),4x,'(',f5.1,'%)',5x,'scf converged')
 1136 format(5x,i4,2(5x,f15.9),4x,'(',f5.1,'%)')
 1140 format(5x,i4,2(5x,f15.9),5x,'scf converged')
 1141 format(5x,i4,2(5x,f15.9))
 1145 format(5x,'one-electron energy',10x,f15.9,
     $      /5x,'coulomb energy     ',10x,f15.9,
     $      /5x,'exchange energy    ',10x,f15.9,
     $      /5x,'correlation energy ',10x,f15.9)
 1190 format (5x,'final vector:')
 1200 format (5x,'punching occupied vectors :')
 1210 format (5(2x,f8.3))
 1215 format (' coulomb matrices ')
 1216 format (' difference coulomb matrices ')
 1220 format (' exchange matrices ')
 1225 format(5x,'exchange-correlation energy',2x,f15.9)
 1226 format(5x,'total integrated charge:',f20.10)
 1230 format(5x,'total becke charges:',/,
     $      (10x,a8,2x,f11.5))
 1240 format(5x,'total becke charges and spins:',/,
     $      (12x,a8,2f11.6))
 1250 format(5x,'time for j-matrix:',18x,f8.1)
 1260 format(5x,'time for k-matrix:',18x,f8.1)

      ierr=stderr()
c
c     --- note t1 and bmatrx are implicitly equivalenced ---
c         set the pointers to diis storage 
      if(dodiis) then
         call izero(dsptr,mxds)
         do 20 i=1,mxiter
            ptr(i)=-9999999
   20    continue
      endif
c
c     --- perhaps using fock matrix damping ---
      if(logkey(ops,'scf=damp',.false.,' ')) then
         call fparr(ops,'scf=damp',damp,nshell-1,' ')
         write(iout,1000) (damp(i),i=1,nshell-1)
      end if
c
c     --- determine options in force ---
      printj=logkey(ops,'scf=print=j',.false.,' ')
      printk=logkey(ops,'scf=print=k',.false.,' ')
      prteexch=logkey(ops,'scf=prteexch',.false.,' ')
c     --- should we abort if convergence not met ---
      abort=.true.
      if(logkey(ops,'scf=noabort',.false.,' ')) then
         abort=.false.
      endif
c
c     --- print t and v if requested ---
      if(logkey(ops,'scf=print=kinetic',.false.,' ')) then
         write (iout,1010)
         call print(t,nnp,nbf,iout)
      end if
      if(logkey(ops,'scf=print=potential',.false.,' ')) then
         write (iout,1020)
         call print(v,nnp,nbf,iout)
      end if
c
c     --- direct scf stuff ---
      if (dirj) then
         cutexp=fpkey(ops,'int=preexponential',1.0d-12,' ')
         qtest=fpkey(ops,'scf=qtest',1.0d-09,' ')
         if(prnt) then
            write(iout,1025) cutexp,qtest
         endif
         rhotest=logkey(ops,'scf=rhotest',.true.,' ')
c        --- make sure we have room and get estimates of integrals.
         call sizer(noprim,nbtype,nocart,maxmom,
     $              natoms,npf,nbf,ndmat,biggy)
         call estmq(ptprim,noprim,nbtype,nocont,ptcont,cont,
     $              ncont,ex,coord,nx,ny,nz,
     $              lenxyz,nocart,mintyp,maxmom,
     $              values,left,natoms,nprim,ops,
     $              valuesi,pstart,nbf,nnp,nnshl,qint)
         if(biggy.gt.left) then
            write(iout,1026)
            write(iout,1027) left,biggy
            call lnkerr('more core for directj') 
         endif
      else if(poisj) then
c
c        --- get the parameters to be used for the atomic decompositions
         vorder=intkey(ops,'poisson=spline-order',3,' ')
         vorder=vorder+1
         vlmax=intkey(ops,'poisson=lmax',4,' ')
         vradial=intkey(ops,'poisson=nradial',51,' ')
         vncrule=intkey(ops,'poisson=rule',5,' ')
         write(iout,1028) vlmax,vradial,vncrule,vorder-1
      endif
c
c     --- form h ---
      call vadd(h,t,v,nnp)
      if(logkey(ops,'scf=print=one-electron',.false.,' ')) then
         write (iout,1030)
         call print(h,nnp,nbf,iout)
      end if
c
c     --- retrieve guess vector and normalize ---
      call iosys('read real "guess vector" from rwf',-1,c,0,' ')
      call iosys('read real "guess orbital energies" from rwf',
     $           -1,eigval,0,' ')
      call schmdt(c,s,t1,t2,t3,nbf,nbf,nnp,maxerr)
      if(logkey(ops,'scf=print=guess',.false.,' ')) then
         write (iout,1040)
         call vecout(c,eigval,nbf)
      end if
      maxolap=logkey(ops,'scf=maxolap',.false.,' ')
      if (maxolap) then
         frcolap=intkey(ops,'scf=frcolap',1,' ')
         write(iout,*)'Starting to force olaps at orbital ',frcolap
      endif
      if (maxolap) then
         call iosys('write real "saved orthog guess" to rwf',nbf*nbf,
     $        c,0,' ')
      endif
c
c     --- get the nuclear repulsion ---
      call iosys('read real "nuclear repulsion energy" from rwf',
     $            1,enuc,0,' ')
      if(prnt) then
         write(iout,1050) enuc
      endif
c
c     --- zero density and fock matrices for last iteration ---
c
c     --- initialize the loop for iterations ---              
      ndiis=0
      iter=0
      mndiis=1
      maxocc=min(nbf,shlmax(nshell-1))
      call rzero(energs,mxiter)
      call rzero(diissv,mxiter)
      if (prnt) then
         if(dirj) then
            write(iout,1065)
         else
            write(iout,1060)
         endif
      end if
      call rzero(dlast,nnp*(nshell-1))
      call rzero(jlast,nnp*(nshell-1))
      call rzero(timej,3)
      call rzero(timek,3)
c
c     --------------------------------------------------------
c     --------------------------------------------------------
c     the iteration loop
      elast=zero
  100 continue
      iter=iter+1
c
c     ----- form the density matrices -----
c     note that the occupation number is not included in the density.
      do 110 i=1,nshell-1
         call gdmat(d(1,i),c,nbf,nnp,shlmin(i),shlmax(i))
  110 continue
      if (debug) then
         write (iout,1070)
         do 120 i=1,ndmat
            call print(d(1,i),nnp,nbf,iout)
  120    continue
      end if
c
c     --- form the j and k matrices ---
      call timing(dum1,dum2,dum3)
      if (dirj) then
c        -- form the difference density matrix.
c           use t1 to hold the difference. note that it is given room
c           equivalent to max(nbf**2,nnp*ncoul) in the calling routine.
c           the density difference in the primitive basis comes back in dpr.
         call vsub(t1,d,dlast,nnp*ndmat)
         do 122 i=1,ndmat
            call tr1dm(values(1),values(1+npf*npf),
     $           t1(1+nnp*(i-1)),prtoao,dpr(1,i),nbf,nnp,npf,
     $           nnprim)
            if(debug) then
               write(iout,*) 'difference primitive density matrices'
               call print(dpr(1,i),nnprim,npf,iout)
            endif
c
c           --- now form an array which contains the largest difference
c               density in the contracted basis for each shell block
c               combination. this will be used in conjunction with the
c               largest estimated integral in the block to screen zeroes.
  122    continue
         call bigdij(noprim,nbtype,nocont,ncont,nocart,natoms,
     $        start,nbf,nnp,t1,nnshl,dijmax,values,left,ndmat)
         call directj(jmat,jprim,ptprim,noprim,nbtype,ex,coord,
     $                nx,ny,nz,lenxyz,nocart,mintyp,maxmom,
     $                values,left,natoms,npf,nnprim,nprim,ops,cutexp,
     $                rhotest,values,.false.,ndmat,
     $                pstart,prtoao,dpr,
     $                nbf,nnp,ntotal(iter),nonzer(iter),nnshl,qint,
     $                dijmax,qtest)
c        --- add the difference j-matrix to the last one ---
         write(iout,1216)
         do 1224 i=1,ncoul
            call print(jmat(1,i),nnp,nbf,iout)
 1224    continue 

         call vadd(jmat,jmat,jlast,nnp*(nshell-1))
      else if (poisj) then
c        call lnkerr('Heh.  nice try')
c        need to set up gofish's core
c        need to set up an itch and an aitch
c        need to allocate and set up total density matrix dtot=2*d for 
c           closed shell
c        dtot=2*d(closed)+d(open) for open shells.
c        need imaxl
c        need to modify gofish to actually do the quadrature too.
c
c        --- compute entire grid for gofish.
         
         dtot=1
         imaxl=wpadti(dtot+nnp)
         itch=iadtwp(imaxl+natoms)
         aitch=wpadti(itch)
         call rzero(values(dtot),nnp)
         call vmove(values(dtot),d,nnp)
c        --- what the hell do we do for open shells?
         pleft=left-nnp-natoms
         if(old) then
            call poisson(values(itch),values(dtot),nbf,nnp,jmat,ncoul,
     $        ndmat,natoms,mxgrd,dmcut,dencut,
     $        mxgbsiz,valuesi(aitch),ian,coord,ex,cont,ptprim,noprim,
     $        nocont,ptcont,mxcont,nprim,ntypes,nbtype,ncont,
     $        start,nocart,nobf,maxmom,minmom,mintyp,nx,ny,nz,
     $        xyzgrid,wts,charge,valuesi(imaxl),bigl,ops,
     $        pleft,rnuc,amu,pwtx,rr,radii,akl,rmax,lmax,nomega,
     $        nradial,grdtyp,adjust,minesz,vorder,vlmax,vradial,
     $        vncrule)
         else
c           some core necessary for gofish2.
            nlm=(vlmax+1)*(vlmax+1)
            rpts=itch
            vlm=rpts+vradial*natoms
            y2=vlm+vradial*nlm*natoms
            itch=y2+vradial*nlm*natoms
            aitch=wpadti(itch)
            pleft=pleft-vradial*natoms*(1+nlm+nlm)
            call rzero(jmat,nnp*ncoul)
c           --- pass the difference density to gofish2.
            call vsub(values(dtot),d(1,1),dlast(1,1),nnp)
            call gofish2(values(itch),values(dtot),nbf,nnp,jmat,ncoul,
     $        ndmat,natoms,mxgrd,dmcut,dencut,
     $        mxgbsiz,valuesi(aitch),ian,coord,ex,cont,ptprim,noprim,
     $        nocont,ptcont,mxcont,nprim,ntypes,nbtype,ncont,
     $        start,nocart,nobf,maxmom,minmom,mintyp,nx,ny,nz,
     $        xyzgrid,wts,charge,valuesi(imaxl),bigl,ops,
     $        pleft,vwts,rnuc,amu,pwtx,rr,radii,akl,ptrad,rmax,lmax,
     $        nomega,nradial,grdtyp,adjust,minesz,vlmax,vradial,
     $        vncrule,nlm,values(rpts),values(vlm),values(y2))
            call vneg(jmat,jmat,nnp)
c           --- and add it to the previous one.
            call vadd(jmat,jmat,jlast,nnp*(nshell-1))
         endif
       
      else
         call jmatrix(values,d,nbf,nnp,jmat,ncoul,ntriang,ndmat)
      endif
      call timing(dum4,dum5,dum6)
      timej(1)=timej(1)+dum4-dum1
      timej(2)=timej(2)+dum5-dum2
      if(timeit) then
         write(iout,*) 'time for j-matrix',
     $                 dum4-dum1,dum5-dum2,dum6-dum3
      endif
c
      if (printj) then
         write(iout,1215)
         do 1222 i=1,ncoul
            call print(jmat(1,i),nnp,nbf,iout)
 1222    continue 
      endif
c
      calce=.true.
      call timing(dum1,dum2,dum3)
      call kmatrix(values,d,dlast,nbf,nnp,kmat,nexch,wts,
     $     ndmat,natoms,mxgrd,exc,slater,becke,lyp,vwn,calce,calc,
     $     dmcut,dencut,eexch,prteexch,
     $     mxgbsiz,valuesi,
     $     coord,ex,cont,ptprim,noprim,nocont,ptcont,
     $     mxcont,nprim,ntypes,nbtype,ncont,
     $     start,nocart,nobf,maxmom,minmom,mintyp,
     $     nx,ny,nz,ian,xyzgrid,charge,maxl,dograd,bigl,
     $     vwts,rnuc,amu,pwtx,rr,radii,akl,ptrad,
     $     rmax,lmax,nomega,nradial,grdtyp,adjust,minesz)
      call timing(dum4,dum5,dum6)
      timek(1)=timek(1)+dum4-dum1
      timek(2)=timek(2)+dum5-dum2
      if(timeit) then
         write(iout,*) ' time for kmatrix',dum4-dum1,dum5-dum2,dum6-dum3
      endif
      if (printk) then
         write(iout,1220)
         do 1221 i=1,nexch
            call print(kmat(1,i),nnp,nbf,iout)
 1221    continue 
         write(iout,1225) exc
      endif

c
c     --- calculate the energy ---
      call ksenrg(jmat,exc,kmat,d,nnp,nbf,ncoul,nexch,ndmat,h,enuc,
     $            energy,calc,eonel,etwoel,
     $            nshell,fcoef,alpha,beta)
      energs(iter)=energy
c
c     --- save the density and coulomb matrix
      call vmove(dlast,d,nnp*(nshell-1))
      call vmove(jlast,jmat,nnp*(nshell-1))
c
c     --- transform the one-electron,j and k matrices to either the
c         pseudo-canonical representation or the current mo basis.
c         if the page-mciver algorithm is in effect, we must transform
c         to the pseudo-canonical representation.
      if(logkey(ops,'scf=pseudo',.false.,' ')
     $  .or.page) then
         call pseudo(nbf,nnp,nshell,ncoul,nexch,h,t1,jmat,kmat,
     $              c,t2,t3,eigval,fcoef,alpha,beta,shlmin,shlmax,t4,t5,
     $              ops)
c        need density matrices in the pseudo-canonical representation
c        in case we are doing diis. 
         do 125 i=1,nshell-1
            call gdmat(d(1,i),c,nbf,nnp,shlmin(i),shlmax(i))
  125    continue
      else
         call tomo(nbf,nnp,nshell,ncoul,nexch,h,t1,jmat,kmat,c,t2,t3)
      endif
c
c     --- form the 'fock' matrices ---
      call qmat(jmat,kmat,t1,nnp,ncoul,nexch,f,pulay,page,
     $          calc,shlmin,shlmax,nshell,nbf,energy,
     $          fcoef,alpha,beta,damp)
      if (logkey(ops,'scf=print=q-matrix',.false.,' ')) then
         write (iout,1013)
 1013    format (/,t10,'q-matrix in the m. o. basis')
         call print(f,nnp,nbf,iout)
      end if
c
c     --- transform the 'fock' matrix back to the ao basis ---
      call trtosq(t3,s,nbf,nnp)
      call trtosq(t1,f,nbf,nnp)
      call ebct(t2,t1,c,nbf,nbf,nbf)
      call ebc(t1,c,t2,nbf,nbf,nbf)
      call ebc(t2,t1,t3,nbf,nbf,nbf)
      call ebc(t1,t3,t2,nbf,nbf,nbf)
      call sqtotr(f,t1,nbf,nnp)
c
c     --- form the total density matrix ---
      call rzero(triang,nnp)
      do 140 shell=1,nshell-1
         do 130 i=1,nnp
            triang(i)=triang(i)+fcoef(shell)*d(i,shell)
  130    continue
  140 continue
c
c     --- perhaps perform pulay's diis extrapolation ---
c         if not, at least compute the error matrix for convergence tests
      if(dodiis) then
         call diis(f,t1,t2,t3,t4,t5,triang,s,smhalf,error,
     $             diissv,bmatrx,focksv,c,diiser,stdiis,
     $             nnp,nbf,mxds,mxiter,iter,ndiis,mndiis,
     $             ptr,dsptr,ops)
      else
         call errmat(f,t1,t2,t3,t4,t5,triang,s,diiser,
     $               smhalf,nnp,nbf,ops)
         diissv(iter)=diiser
      endif
c
c     --- transform back to the mo basis, yet again ---
      call trtosq(t2,f,nbf,nnp)
      if (usesym) then
ctemp        call ebc(t1,t2,salc,nbf,nbf,nbf)
ctemp        call ebtc(t2,salc,t1,nbf,nbf,nbf)
      else
         call ebc(t1,t2,c,nbf,nbf,nbf)
         call ebtc(t2,c,t1,nbf,nbf,nbf)
      endif
      call sqtotr(f,t2,nbf,nnp)
c
c     --- perhaps apply a level-shift to the virtual orbitals ---
      if (level.ne.zero) then
         call levshft(f,nnp,nbf,maxocc,level)
      endif
c
c     --- print the fock matrix ---
      if (debug) then
         write (iout,1090)
         call print(f,nnp,nbf,iout)
      end if
c
c     --- move a copy of the fock matrix into triang ---
      call vmove(triang,f,nnp)
c
c     --- diagonalize f ---
      if(usesym) then
         if (logkey(ops,'scf=print=salc-fock-mtx',.false.,' '))then
            call trtosq(t2,triang,nbf,nnp)
            call ebc(t3,t2,salc,nbf,nbf,nbf)
            call ebtc(t3,salc,t3,nbf,nbf,nbf)
            call sqtotr(t5,t3,nbf,nnp)
            write(iout,*)'The fock matrix in the salc basis'
            call print(t5,nnp,nbf,iout)
         endif
         call symrsp(triang,u,t2,t3,t4,t5,eigval,temp,kept,nbf,nnp,
     $               nshell,shlmin,shlmax,nirrep,numso,lambda,
     $               orbsym,occsym,salc,c)
c        fill in the orbital symmetry labels.
         do 150 i=1,nbf
            symlabl(i)=lirrep(orbsym(i))
  150    continue
      else
         call degrsp(nbf,nnp,triang,eigval,1,u,t2,t3)
c        --- form new scf vector = c * u ---
         call ebc(t1,c,u,nbf,nbf,nbf)
         call vmove(c,t1,nbf**2)
      endif
c
c     --- determine orbital occupancy ---
      if (maxolap)then
c        occupy according to maximal overlap with the guess vector.
c        form guess^T S c         
         call trtosq(t2,s,nbf,nnp)
         call ebc(t1,t2,c,nbf,nbf,nbf)
         call iosys('read real "saved orthog guess" from rwf',nbf*nbf,
     $               t3,0,' ')
         call ebtc(t2,t3,t1,nbf,nbf,nbf)
         if (logkey(ops,'scf=maxolap=debug',.false.,' ')) then
            write(iout,*)'Overlap between current and guess'
            call vecout(t2,eigval,nbf)
         endif
c
c        now try to flip stuff around until we have roughly 1 on the diags.
c        loop over guess orbitals, check until find largest
         do 155 i=frcolap,maxocc
            tmpmx=0.0
            tmpindx=0
            do 154 j=1,nbf
               if (abs(t2(i,j)).gt.tmpmx) then
                  tmpmx=abs(t2(i,j))
                  tmpindx=j
               endif
 154        continue 
            if (tmpindx.eq.i) goto 155
            if (tmpindx .eq.0) call lnkerr(
     $           'ack: cannot find large olap')
c
c              swap the i'th new orbital with the tmpindx'th
            call vmove(t1,c(1,i),nbf)
            call vmove(c(1,i),c(1,tmpindx),nbf)
            call vmove(c(1,tmpindx),t1,nbf)
            call vmove(t1,t2(1,i),nbf)
            call vmove(t2(1,i),t2(1,tmpindx),nbf)
            call vmove(t2(1,tmpindx),t1,nbf)
            tmpmx=eigval(i)
            eigval(i)=eigval(tmpindx)
            eigval(tmpindx)=tmpmx
            if (usesym) then
               itmp=orbsym(i)
               orbsym(i)=orbsym(tmpindx)
               orbsym(tmpindx)=itmp
               symlabl(i)=lirrep(orbsym(i))
               symlabl(tmpindx)=lirrep(orbsym(tmpindx))
            endif
 155     continue 
         if (logkey(ops,'scf=maxolap=debug',.false.,' ')) then
            call trtosq(t2,s,nbf,nnp)
            call ebc(t1,t2,c,nbf,nbf,nbf)
            call ebtc(t2,t3,t1,nbf,nbf,nbf)
            write(iout,*)'After swapping'
            call vecout(t2,eigval,nbf)
         endif
      endif
c
c     --- print rotation matrix ---
      if (debug) then
         write (iout,1120)
         call vecout(u,eigval,nbf)
      end if
c
c     --- print vectors ---
      if (logkey(ops,'scf=print=viters',.false.,' ')) then
         write (iout,1130) iter
         if(usesym) then
            write(iout,*) ' orbital symmetries'
            write(iout,*) (orbsym(i),i=1,nbf)
         end if
         call vecout(c,eigval,nbf)
      end if
c
c     --- convergence test ---
c         first check the diis convergence.
      if (diiser.lt.ten**(-finish)) go to 200
c         now check change in energy. 
c         with the smaller integration grids it is not always possible
c         to reach the default diis accuracies.
c         make sure that the energy change is an order of magnitude smaller.
      deltae=abs(energs(iter)-elast)
      elast=energs(iter)
      if(deltae.le.ten**(-finish-1)) go to 200
      if(dirj) then
         write(ierr,1136) iter,energs(iter),diissv(iter),
     $            float(nonzer(iter))*100.d0/float(ntotal(iter))
      else
         write(ierr,1141) iter,energs(iter),diissv(iter)
      endif

      if (iter.lt.mxiter) go to 100
c     --- end iteration loop ---
c     -------------------------------------------------------------------
c
c     --- convergence not met in mxiter iterations ---
      do 190 i=1,iter
         if(dirj) then
            write(iout,1136) i,energs(i),diissv(i),
     $            float(nonzer(i))*100.d0/float(ntotal(i))
         else
            write(iout,1141) i,energs(i),diissv(i)
         endif
  190 continue
      if(abort) then
         call lnkerr('convergence not met.')
      else
         call iosys('write real "scf unconverged vector" on rwf',
     $               nbf**2,c,0,' ')
         call iosys('write real "scf unconverged orbital energies"'
     $               //' on rwf',nbf,eigval,0,' ')
      end if
c
  200 continue
      if (debug) then
         do 290 i=1,iter
            if(dirj) then
               write(iout,1136) i,energs(i),diissv(i),
     $            float(nonzer(i))*100.d0/float(ntotal(i))
           else
               write(iout,1141) i,energs(i),diissv(i)
            endif
  290    continue
      endif
c
c     --- print the final energy and vector if desired ---
      if (prnt.and..not.logkey(ops,'print=scf=energy',.false.,' ')) then
         if(dirj) then
            write(iout,1135) iter,energs(iter),diissv(iter),
     $            float(nonzer(iter))*100.d0/float(ntotal(iter))
         else
            write(iout,1140) iter,energs(iter),diissv(iter)
         endif
      endif
      if (prteexch) then
         write(iout,1145) eonel,etwoel-exc,eexch,exc-eexch
      endif
      if (logkey(ops,'print=m511=vector',.false.,' ')) then
         call iosys('read character "basis function labels" from rwf',
     $               len(bflabl(1))*nbf,0,0,bflabl)
         write (iout,1190)
         if (chrkey(ops,'print=m511=vector',' ',' ').eq.'all') then
            nprint=nbf
         else
            nprint=min(maxocc+4,nbf)
         end if
         if(usesym) then
            call matprt(c,nbf,nbf,nbf,nprint,1,2,bflabl,symlabl,
     $                  0,eigval,.true.)
         else
            call wvec(c,eigval,nbf,nprint,bflabl,' ')
         endif
      end if
c
c     --- determine the total integrated charge.
      totchg=zero
      do 300 i=1,natoms
         if(calc.eq.'closed') then
            totchg=totchg+charge(i,1)
         else if(calc.eq.'open') then
            totchg=totchg+charge(i,1)+charge(i,2)
         endif
  300 continue
      if(prnt) then
         write(iout,1226) totchg
      endif
c
c     --- print the atomic charges according to the becke
c         fuzzy voronoi polyhedra scheme if desired
      if(prnt.and.logkey(ops,'print=scf=charges',.false.,' ')) then
         if(calc.eq.'closed') then
            write(iout,1230) (names(i),charge(i,1),i=1,natoms)
         else if(calc.eq.'open') then
            write(iout,1240) (names(i),charge(i,1)+charge(i,2),
     $                        charge(i,1)-charge(i,2),i=1,natoms)
         endif
      endif
c
c     --- print the vectors in a manner in which m401 can read them
c     from cards with the rdinp option.
      if(logkey(ops,'scf=punch',.false.,' '))then
         write(iout,1200)
         write(iout,*)'$vectors'
         do 210 i=1,maxocc
            write(iout,*)' vector ',i
            write(iout,1210) (c(j,i),j=1,nbf)
 210     continue
         write(iout,*)'$end'
      end if
c
c     --- store the energy, vector, eigenvalues, and density
c         matrices
      call iosys('write real f to rwf',nshell,fcoef,0,' ')
      call iosys('write real energy to rwf',1,energy,0,' ')
      call iosys('write real "hf energy" to rwf',1,energy,0,' ')
      call iosys('write real "hf 1e energy" to rwf',1,eonel,0,' ')
      call iosys('write real "hf 2e energy" to rwf',1,etwoel,0,' ')
      ecoul=etwoel-exc
      call iosys('write real "hf coulomb energy" to rwf',
     $            1,ecoul,0,' ')
      call iosys('write real "hf xc energy" to rwf',
     $            1,exc,0,' ')
      call putvec(c,eigval,nbf,usesym,orbsym,symlabl)
      do 220 i=1,nshell-1
         call gdmat(d(1,i),c,nbf,nnp,shlmin(i),shlmax(i))
 220  continue
      call iosys('write real "hf density matrix" to rwf',nnp*ndmat,
     $            d,0,' ')
c
c     --- report total times
      if(prnt) then
         write(iout,1250) timej(1)
         write(iout,1260) timek(1)
      endif
c
c
      return
      end