*deck %W%  %G%
      subroutine dxopen(dengrida,dengridb,dengrada,dengradb,ga,gb,gab,
     $     fout,itch,phi,grad,scr,scr2,slater,becke,lyp,vwn,d,
     $     nbf,nnp,kmat,nexch,grdwt,ndmat,nat,mxgrd,calce,exc,ngrid,
     $     dmcut,dencut,calc,nzptrs,tmpgwt,eexch,
     $     prteexch,queue,tea,kay,ktmp,dorhoup,ngb,gblksiz,mxgbsiz,
     $     mxgblk,c,ex,cont,ptprim,noprim,nocont,ptcont,
     $     mxcont,nprim,ntypes,nbtype,ncont,
     $     start,nocart,nobf,maxmom,minmom,mintyp,
     $     nx,ny,nz,xyzgrid,charge,maxl,dograd,bigl)
c***begin prologue     dxopen.f
c***date written       yymmdd  
c***revision date      2/5/95      
c   oct 24, 1995       russo at lanl
c      finally making open shell gradients actually work
c   may 13, 1994       rlm at lanl
c      modifying kmopen(m511) to do derivatives of exchange-correlation
c      energy.
c***keywords           
c***author             russo, thomas(lanl)
c***source             %W%   %G%
c***purpose            
c***description
c     
c    
c
c***references
c
c***routines called
c
c***end prologue       dxopen.f
      implicit none
c     --- input variables -----
c     --- input arrays (unmodified) ---
c     --- input arrays (scratch) ---
c     --- output arrays ---
c     --- output variables ---
c     --- scratch arrays ---
c     --- local variables ---
      integer nat,mxgrd,ndmat,nexch,nnp,nbf,mxgbsiz,mxgblk
      integer nprim,ntypes,nbtype,ncont,mxcont,bigl
      integer ngrid(nat),ngb(nat),gblksiz(mxgblk,nat)
      real*8 slalph,beckb,lypa,lypb,lypc,lypd,dmcut,dencut
      real*8 dengrida(mxgbsiz),dengridb(mxgbsiz),dengrada(mxgbsiz,3),
     $     dengradb(mxgbsiz,3),ga(mxgbsiz),gb(mxgbsiz),gab(mxgbsiz),
     $     fout(mxgbsiz,*),itch(*),phi(mxgbsiz,nbf),grad(mxgbsiz,nbf,3),
     $     scr(mxgbsiz),scr2(mxgbsiz)
      logical slater,becke,lyp,vwn,calce,prteexch,dorhoup
      logical directk,dograd
      real*8 d(nnp,ndmat),kmat(nnp,nexch),exc,grdwt(mxgrd,nat),eexch
      real*8 charge(nat,ndmat)
      integer ptprim(nat,ntypes),noprim(nat,ntypes),nocont(nat,ntypes)
      integer ptcont(nat,ntypes),start(nat,ntypes),nocart(0:*)
      integer nobf(ntypes),maxmom(ntypes),minmom(ntypes),mintyp(ntypes)
      integer nx(*),ny(*),nz(*)
      real*8 c(3,nat),ex(nprim),cont(ncont)
      real*8 xyzgrid(mxgrd,3,nat)
      integer x,t,omega,delta,iatom,mu,nu,kl,derivs,coord
      integer rho,zeta,g,h,t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11
      real*8 zero,half,one,two,four
      real*8 sdot,foo1,foo2
      real*8 tmpgwt(mxgbsiz)
      real*8 queue(mxgbsiz,3),tea(mxgbsiz,nbf),kay(nbf,nbf),ktmp(nnp)
      integer maxl(nat)
      integer nzptrs(mxgbsiz)
      integer ngtmp,i,bsetptr,bsgptr,gradptr,denptr,ioff,ng,igblk
      integer xyzpow,r,s,rsq,tbf,hess
      real*8 dum1,dum2,dum3,dum4,dum5,dum6
      character*(*) calc
      integer inp,iout
      integer top,ngot,wpadti
      logical timeit
      real*8 timbf,timden,timsqz,timfcn,timk
c
      data timbf,timden,timsqz,timfcn,timk/5*0.0d+00/
      save timbf,timden,timsqz,timfcn,timk
c
      integer hescod(3,3)
      data hescod/1,4,5,4,2,6,5,6,3/

      common /io/ inp,iout
c
      parameter (timeit=.false.)
      parameter (zero=0.0d+00,half=0.5d+00,one=1.0d+00,two=2.0d+00)
      parameter (four=4.0d0)
      parameter (slalph=2.0d0/3.0d0,beckb=0.0042d0,lypa=0.04918d0)
      parameter (lypb=0.132d0,lypc=0.2533d0,lypd=0.349d0)
c
 1010 format(1x,80a1)
c
      call getscm(0,itch,ngot,'dxopen',0)
c
c     --- core allocation within itch
c         note that the same region is used for scratch in both
c         bfgrd and the functional routines.
c    
c     --- loop over atoms, calculate gradient contributions from each
c
      next=nxtval(nproc)+1
      do 10 iatom=1,nat
         if (iatom .ne. next) goto 10
c
c         ---make the grid
c
         call mkatmg(c,ian,xyzgrid,grdwt,rmax,lmax,nomega,nradial,
     $        nat,ngrid,mxgrd,vwts,rnuc,amu,pwtx,rr,adjust,radii,
     $        akl,grdtyp(iatom),radshls,ptrad,.true.,gradwts,iatom)
c
c         ---partition into grid blocks
c
         ngb=ngrid/mxgbsiz
         oddblock=mod(ngrid,mxgbsiz)
         if (oddblock.ne.0) ngb=ngb+1
         ioff=0
         call rzero(dkay,nbf*nbf*3*2)
c
c         ---loop over grid blocks
c
         do 20 igblk=1,ngb
            if (igblk .ne. ngb .or. oddblock.eq.0) then
               ng=mxgbsiz
            else
               ng=oddblock
            endif
c
c            ---get basis, gradients and hessians on this grid block
c
            call bfgrd(c,ex,cont,ptprim,noprim,nocont,ptcont,
     $           nat,nprim,ntypes,nbtype,nnp,ncont,
     $           start,nbf,nocart,nobf,maxmom,minmom,mintyp,
     $           nx,ny,nz,xyzgrid,mxgrd,
     $           itch(xyzpow),scr,itch(rsq),itch(s),itch(r),
     $           itch(tee),phi,grad,hess,.true.,.true.,mxcont,mxgbsiz,
     $           ng,ioff,maxl)
c
c            ---get density and its gradient on this block
c
            call rzero(dengrida,mxgbsiz)
            call rzero(dengridb,mxgbsiz)
            if (fgrad) then
               call rzero(dengrada,mxgbsiz*3)
               call rzero(dengradb,mxgbsiz*3)
            endif
            call gridden(nbf,nnp,ng,mxgbsiz,d(1,1),phi,grad,
     $           dengridb,dengradb,minesiz,scr,itch(phibar),
     $           itch(gradbar),dmcut,fgrad)
            call gridden(nbf,nnp,ng,mxgbsiz,d(1,2),phi,grad,
     $           dengrida,dengradb,minesiz,scr,itch(phibar),
     $           itch(gradbar),dmcut,fgrad)
            call vadd(dengrida,dengrida,dengridb,ng)
            if (fgrad) then
               call vadd(dengrada(1,1),dengrada(1,1),dengradb(1,1),
     $              ng)
               call vadd(dengrada(1,2),dengrada(1,2),dengradb(1,2),
     $              ng)
               call vadd(dengrada(1,3),dengrada(1,3),dengradb(1,3),
     $              ng)
            endif
c
c            ---squeeze based on total density
c
            call sqzopn(ng,ngtmp,mxgbsiz,nbf,phi,grad,dengrida,
     $           dengrada,dengridb,dengradb,grdwt(1+ioff),
     $           tmpgwt,nzptrs,dencut,becke,lyp)
c
c           ---form gradient invariants
c
            if (fgrad) then
               call vmul(gb,dengradb,dengradb,ngtmp)
               call vwxy(gb,gb,dengradb(1,2),dengradb(1,2),+1,
     $              ngtmp)
               call vwxy(gb,gb,dengradb(1,3),dengradb(1,3),+1,
     $              ngtmp)
               call vmul(ga,dengrada,dengrada,ngtmp)
               call vwxy(ga,ga,dengrada(1,2),dengrada(1,2),+1,
     $              ngtmp)
               call vwxy(ga,ga,dengrada(1,3),dengrada(1,3),+1,
     $              ngtmp)
               if (lyp) then
                  call vmul(gab,dengrada,dengradb,ngtmp)
                  call vwxy(gab,gab,dengrada(1,2),dengradb(1,2),+1,
     $                 ngtmp)
                  call vwxy(gab,gab,dengrada(1,3),dengradb(1,3),+1,
     $                 ngtmp)
                  call vclean(gab,dencut,ngtmp)
               endif
               call vclean(ga,dencut,ngtmp)
               call vclean(gb,dencut,ngtmp)
            endif
c
c            ---get functionals/derivatives
c
c remember: fout(.,1)=f, fout(.,2)=df/drhoa, fout(.,3)=df/dgammaaa
c                        fout(.,4)=df/drhob, fout(.,5)=df/dgammabb
c              fout(.,6)=df/dgammaab
c
            derivs=1
            if (slater) then
               call slaterf(ngtmp,mxgbsiz,dengrida,dengridb,slalph,
     $              derivs,fout)
            else if (becke) then
c              still have that silly 'uhf' garbage in there.
               call beckef(ngtmp,mxgbsiz,dengrida,dengridb,ga,gb,
     $              itch(x),itch(t),derivs,beckb,fout,'uhf',1)
            endif
            if (prteexch) then
               eexch=eexch+sdot(ngtmp,tmpgwt,1,fout,1)
            endif
            if (vwn) then
               call vwnf(ngtmp,mxgbsiz,dengrida,dengridb,calc,derivs,
     $              fout,itch(rho),itch(x),itch(zeta),itch(g),
     $              itch(h),itch(t1),itch(t2),itch(t3),itch(t4),
     $              itch(t5),itch(t6),itch(t7),itch(t8),itch(t9),
     $              itch(t10),itch(t11))
            else if (lyp) then
               call rzero(fout(1,6),mxgbsiz)
               call lypf(ngtmp,mxgbsiz,dengrida,dengridb,ga,gb,gab,
     $              itch(omega),itch(delta),itch(t),derivs,lypa,
     $              lypb,lypc,lypd,fout)
            endif
            exc=exc+sdot(ngtmp,tmpgwt,1,fout,1)
c
c            ---do weight-derivative calculation
c
            do 6060 jatom=1,nat
               do 6061 coord=1,3
                  call gthr(scr,gradwts(ioff+1,coord,jatom),
     $                 nzptrs,ngtmp)
                  g1=sdot(ngtmp,fout,1,scr,1)
                  grade(coord,jatom)=grade(coord,jatom)+g1
                  grade(coord,iatom)=grade(coord,iatom)-g1
 6061          continue 
 6060       continue 
c
c            ---form terms independent of gradient-corrections
c
            call vmul(scr,fout(1,2),tmpgwt,ngtmp)
            do 60 nu=1,nbf
               call vmul(tea(1,nu),scr,phi(1,nu),ngtmp)
 60         continue 
            do 70 coord=1,3
               call sgemm('t','n',nbf,nbf,ngtmp,one,grad(1,1,coord),
     $              mxgbsiz,tea,mxgbsiz,one,dkay(1,1,coord,1),nbf)
 70         continue 
            call vmul(scr,fout(1,4),tmpgwt,ngtmp)
            do 601 nu=1,nbf
               call vmul(tea(1,nu),scr,phi(1,nu),ngtmp)
 601        continue 
            do 701 coord=1,3
               call sgemm('t','n',nbf,nbf,ngtmp,one,grad(1,1,coord),
     $              mxgbsiz,tea,mxgbsiz,one,dkay(1,1,coord,2),nbf)
 701        continue 
c
c            ---form gradient-corrected terms for functionals that need it
c
c   X(mu,nu)=phi(nu)Grad(Grad phi(mu))^T+Grad(phi(mu))Grad(phi(nu))^T
c         This is a matrix of *MATRICES*  (Grad(Grad phi(mu))^T is a hessian)
c
c  The terms we need to form are then
c
c     dk(mu,nu,1)=sum(grid) [w*X(mu,nu)*(2f3*grad(rhoa)+f6*grad(rhob))]
c and 
c     dk(mu,nu,2)=sum(grid) [w*X(mu,nu)*(f6*grad(rhoa)+2f5*grad(rhob))]

 20      continue 
c     
c         --- do unrestricted nu summation
c
c         --- do restricted sums and put terms into energy gradient.
         next=nxtval(nproc)+1
 10   continue 
c do a barrier before returning
      next=nxtval(-nproc)
      return
      end
