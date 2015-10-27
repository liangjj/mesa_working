*deck @(#)kmtxguts.f	5.1  11/28/95
      subroutine kmtxguts(dengrida,dengrada,ga,
     $     fout,itch,phi,grad,scr,scr2,slater,becke,lyp,vwn,hf,d,nbf,
     $     nnp,
     $     kmat,nexch,grdwt,ndmat,nat,mxgrd,calce,exc,
     $     dmcut,dencut,calc,nzptrs,tmpgwt,eexch,prteexch,queue,tea,
     $     dorhoup,mxgbsiz,
     $     c,ex,cont,ptprim,noprim,nocont,ptcont,
     $     mxcont,nprim,ntypes,nbtype,ncont,
     $     start,nocart,nobf,maxmom,minmom,mintyp,
     $     nx,ny,nz,ian,xyzgrid,charge,maxl,dograd,bigl,
     $     vwts,rnuc,amu,pwtx,rr,radii,akl,ptrad,
     $     rmax,lmax,nomega,nradial,grdtyp,adjust,minesz,phibar,
     $     gradbar,nocc,coeff,bfcut,kmcut,grdfil)
c
c***begin prologue     kmtxguts.f
c***date written       930521   (yymmdd)  
c***revision date      11/28/95      
c
c***keywords           k matrix coulomb matrix closed shell
c***author             RUSSO, thomas (lanl)
c***source             @(#)kmtxguts.f	5.1   11/28/95
c***purpose            to form the k matrix given basis and gradients
c                      on grid, and density matx for the closed shell 
c                      special case
c***description
c     
c    
c
c***references
c
c***routines called
c
c***end prologue       kmtxguts.f
      implicit none
c     --- input variables -----
      integer nat,mxgrd,ndmat,nexch,nnp,nbf,mxgbsiz
      integer nprim,ntypes,nbtype,ncont,mxcont
      integer rmax,lmax,nomega,nradial
      integer minesz,nocc
      real*8 dencut,dmcut
      logical slater,becke,lyp,vwn,calce,prteexch,dorhoup,hf
      logical dograd,adjust
      character*(*) calc
      character*4 grdfil
c     --- input arrays (unmodified) ---
      integer ian(nat)
      integer ptprim(nat,ntypes),noprim(nat,ntypes),nocont(nat,ntypes)
      integer ptcont(nat,ntypes),start(nat,ntypes),nocart(0:*)
      integer nobf(ntypes),maxmom(ntypes),minmom(ntypes),mintyp(ntypes)
      integer nx(*),ny(*),nz(*)
      real*8 c(3,nat),ex(nprim),cont(ncont)
      real*8 d(nnp,nexch)
      real*8 coeff(nbf,nbf)
      character*(*) grdtyp(nat)
c     --- input arrays (scratch) ---
      integer nzptrs(mxgbsiz)
      integer maxl(nat)
      integer ptrad(rmax)
      real*8 dengrida(mxgbsiz),dengrada(mxgbsiz,3),
     $       ga(mxgbsiz),
     $       fout(mxgbsiz,*),itch(*),phi(mxgbsiz,nbf),
     $       grad(mxgbsiz,nbf,3),scr(mxgbsiz),scr2(mxgbsiz)
      real*8 phibar(nbf),gradbar(nbf)
      real*8 grdwt(mxgrd)
      real*8 tmpgwt(mxgbsiz)
      real*8 queue(minesz,3),tea(minesz,nbf)
      real*8 vwts(mxgrd),rnuc(nat,nat),amu(nat,nat)
      real*8 pwtx(nat),rr(nat),radii(nat),akl(nat,nat)
      real*8 xyzgrid(mxgrd,3)
c     --- output arrays ---
      real*8 kmat(nnp,nexch)
      real*8 charge(nat,ndmat)
c     --- output variables ---
      real*8 exc,eexch
c     --- scratch arrays ---
c     --- local variables ---
      integer x,t,omega,delta,iatom,derivs
      integer ngtmp,igblk,ioff,ng
      integer rho,zeta,g,h,t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11
      integer gradptr,denptr
      integer inp,iout
      integer xyzpow,r,s,tbf,hess,rsq
      integer top,ngot,wpadti
      integer bigl
      integer radshls
      integer oddblock,ngb,ngrid
      integer i
      logical dogrdio
      logical new
      real*8 zero,half,one,two,four
      real*8 sdot
      real*8 slalph,beckb,lypa,lypb,lypc,lypd
      real*8 dum1,dum2,dum3,dum4,dum5,dum6
      real*8 timbf,timden,timsqz,timfcn,timk
      real*8 bfcut,kmcut
      logical timeit
c
      data timbf,timden,timsqz,timfcn,timk/5*0.0d+00/
      save timbf,timden,timsqz,timfcn,timk
c
      common /io/ inp,iout
c
      parameter (timeit=.true.)
      parameter (dogrdio=.true.)
      parameter (new=.true.)
      parameter (zero=0.0d+00,half=0.5d+00,one=1.0d+00,two=2.0d+00)
      parameter (four=4.0d0)
      parameter (slalph=2.0d0/3.0d0,beckb=0.0042d0,lypa=0.04918d0)
      parameter (lypb=0.132d0,lypc=0.2533d0,lypd=0.349d0)
c
c TCGMSG
c
      integer mynodeid,nodeid,nproc,nnodes,nxtval,next
      include 'msgtypesf.h'
c
 1010 format(5x,80a1)
c
c
c 
      mynodeid=nodeid()
      nproc=nnodes()
c
c if only HF, and no corrf, this routine does nothing!
c will have to work out hybrid stuff, though.
c
      if (hf .and. (.not. (vwn .or. lyp))) return
c
c     --- core allocation within itch
c         note that the space is used by both bfgrd and the functional
c         routines
c         make sure we can really fit in as much as we think
      call getscm(0,itch,ngot,'kmtxguts',0)
      s=1
      r=s+mxgbsiz*mxcont
      tbf=r
      hess=tbf
      rsq=r+mxgbsiz*mxcont
      xyzpow=rsq+mxgbsiz
      top=xyzpow+3*mxgbsiz*bigl
      if(wpadti(top).gt.ngot) then
         if (mynodeid.eq.0) write(iout,*) 'need,have',wpadti(top),ngot
         call plnkerr('need more memory in kmtxguts',1000)
      endif
c
c     --- functional scratch space
      t=0
      if (slater) then
c        needs nothing
         top=1
      else if (becke) then
         x=1
         t=x+mxgbsiz
         top=t+4*mxgbsiz
      else if (hf) then
c does nothing!
         top=1
      endif
      if (vwn) then
         x=1
         rho=x+mxgbsiz
         zeta=rho+mxgbsiz
         g=zeta+mxgbsiz
         h=g+mxgbsiz
         t1=h+mxgbsiz
         t2=t1+mxgbsiz
         t3=t2+mxgbsiz
         t4=t3+mxgbsiz
         t5=t4+mxgbsiz
         t6=t5+mxgbsiz
         t7=t6+mxgbsiz
         t8=t7+mxgbsiz
         t9=t8+mxgbsiz
         t10=t9+mxgbsiz
         t11=t10+mxgbsiz
         top=t11+mxgbsiz
      else if (lyp) then
         if (t.eq.0) t=1
         omega=t+mxgbsiz*10
         delta=omega+mxgbsiz
         top=top+mxgbsiz
      endif
      if(wpadti(top).gt.ngot) then
         if (mynodeid.eq.0) write(iout,*) 'need,have',wpadti(top),ngot
         call plnkerr('need more core in kmtxguts',1001)
      endif
c     
c     --- loop over atomic grids
      if(dogrdio) then
c        can use the scratch area ptrad to store grid sizes
         call iosys('read "external atomic grid size" from '
     $        //grdfil(1:3),nat,ptrad,0,' ') 
      endif
      call rzero(charge,ndmat*nat)
      denptr=0
      gradptr=0
      next=nxtval(nproc)+1
      do 10 iatom=1,nat
         if (iatom .ne. next) goto 10
c         write(6,*)" Node:",mynodeid," doing atom ",iatom
c
c        --- generate grid, calc number of grid blocks
         if(dogrdio) then
            ngrid=ptrad(iatom)
            call iosys('read real "external grid" from '//grdfil(1:3),
     $                 mxgrd*3,xyzgrid,mxgrd*3*(iatom-1),' ')
            call iosys('read real "external grid weights" from '
     $           //grdfil(1:3),mxgrd,grdwt,mxgrd*(iatom-1),' ')
         else
            call mkatmg(c,ian,xyzgrid,grdwt,rmax,lmax,nomega,
     $                  nradial,
     $                  nat,ngrid,mxgrd,vwts,rnuc,amu,pwtx,rr,adjust,
     $                  radii,akl,grdtyp(iatom),
     $                  radshls,ptrad,.false.,vwts,iatom)
         endif
         ngb=ngrid/mxgbsiz
         oddblock=mod(ngrid,mxgbsiz)
         if (oddblock.ne.0) ngb=ngb+1
         ioff=0
         do 15 igblk=1,ngb
            if (igblk .ne. ngb .or. oddblock.eq.0 ) then
               ng=mxgbsiz
            else
               ng=oddblock
            endif
c
c           --- compute basis functions on this grid block.
            if(timeit) then
              call timing(dum1,dum2,dum3)
            endif
            call bfgrd(c,ex,cont,ptprim,noprim,nocont,ptcont,
     $                 nat,nprim,ntypes,nbtype,nnp,ncont,
     $                 start,nbf,nocart,nobf,maxmom,minmom,mintyp,
     $                 nx,ny,nz,xyzgrid,mxgrd,
     $                 itch(xyzpow),scr,itch(rsq),itch(s),itch(r),
     $                 itch(tbf),phi,grad,itch(hess),dograd,.false.,
     $                 mxcont,mxgbsiz,ng,ioff,maxl,bfcut)
            if(timeit) then
               call timing(dum4,dum5,dum6)
               timbf=timbf+dum4-dum1
            endif
c
c           --- if d contains the full density matrix, 
c               zero dengrid and dengrad.
c               if d contains the difference between this iteration's
c               and the last, load up the last iteration's density.
            if (.not. dorhoup) then
               call rzero(dengrida,mxgbsiz)
               if (becke .or. lyp) then
                  call rzero(dengrada,mxgbsiz*3)
               endif
            else
               call plnkerr('Yo!  Who changed dorhoup?',1003)
               call iosys('read real "closed density" from rwf',
     $                     mxgbsiz,dengrida,denptr,' ')
               if (becke .or. lyp) then
                  call iosys('read real "closed density gradient"'
     $                     //' from rwf',mxgbsiz*3,dengrada,gradptr,' ')
               endif
            endif
c
c           --- form density and gradients on the grid for this atom. 
            if(timeit) then
               call timing(dum1,dum2,dum3)
            endif
            if (new) then
               call moden(nbf,nocc,ng,mxgbsiz,coeff,phi,grad,
     $              dengrida,dengrada,minesz,scr,scr2,
     $              phibar,gradbar,bfcut,dograd)
            else
               call gridden(nbf,nnp,ng,mxgbsiz,d(1,1),phi,grad,
     $                   dengrida,dengrada,minesz,scr,phibar,gradbar,
     $                   dmcut,dograd)
            endif
            if(timeit) then
               call timing(dum4,dum5,dum6)
               timden=timden+dum4-dum1
            endif
c
c$$$            call iosys('write real "closed density" on rwf',
c$$$     $                  mxgbsiz,dengrida,denptr,' ')
c$$$            if (becke .or. lyp) then
c$$$               call iosys('write real "closed density gradient" on rwf',
c$$$     $                     mxgbsiz*3,dengrada,gradptr,' ')
c$$$            endif
            denptr=denptr+mxgbsiz 
            gradptr=gradptr+mxgbsiz*3
c
c
            if(timeit) then
               call timing(dum1,dum2,dum3)
            endif
c
c           --- now we clean the density, and then use nzident 
c               to squeeze out all those values which are zero.
c               given the non-zero density points, squeeze the
c               gradient of the density, the functions and their
c               gradients, and the weights.
            call sqzclo(ng,ngtmp,mxgbsiz,nbf,phi,grad,dengrida,
     $                  dengrada,grdwt(1+ioff),tmpgwt,nzptrs,
     $                  dencut,becke,lyp) 
c
c           --- do no more if ngtmp=0
            if (ngtmp.eq.0) goto 14
c
c           --- form density invariants
            if(becke.or.lyp) then
               call vmul(ga,dengrada,dengrada,ngtmp)
               call vwxy(ga,ga,dengrada(1,2),dengrada(1,2),+1,
     $                   ngtmp)
               call vwxy(ga,ga,dengrada(1,3),dengrada(1,3),+1,
     $                   ngtmp)
               call vclean(ga,dencut,ngtmp)
            endif
            if(timeit) then
               call timing(dum4,dum5,dum6) 
               timsqz=timsqz+dum4-dum1
            endif
c
c           --- sum the density to obtain atomic charges ---
            charge(iatom,1)=charge(iatom,1)
     $                     +two*sdot(ngtmp,dengrida,1,tmpgwt,1)
c
c           --- form functional's derivatives on the grid using that 
c               stuff. both of the exchange functionals clobber the
c               contents of fout.
            if(timeit) then
               call timing(dum1,dum2,dum3)
            endif
            derivs=-1
            if (calce) derivs=1
            if (lyp) call rzero(fout(1,6),mxgbsiz)
            if (slater) then
               call slaterf(ngtmp,mxgbsiz,dengrida,dengrida,slalph,
     $                      derivs,fout)
            else
               call beckef(ngtmp,mxgbsiz,dengrida,dengrida,ga,ga,
     $                     itch(x),itch(t),derivs,beckb,fout,calc,1)
            endif
            if (calce .and. prteexch) then
               eexch=eexch+sdot(ngtmp,tmpgwt,1,fout,1)
            endif
            if (vwn) then
               call vwnf(ngtmp,mxgbsiz,dengrida,dengrida,calc,derivs,
     $                   fout,itch(rho),itch(x),itch(zeta),itch(g),
     $                   itch(h),itch(t1),itch(t2),itch(t3),itch(t4),
     $                   itch(t5),itch(t6),itch(t7),itch(t8),itch(t9),
     $                   itch(t10),itch(t11))
            else if (lyp) then
               call lypf(ngtmp,mxgbsiz,dengrida,dengrida,ga,ga,ga,
     $                   itch(omega),itch(delta),itch(t),derivs,lypa,
     $                   lypb,lypc,lypd,fout)
            endif
            if(timeit) then
               call timing(dum4,dum5,dum6)
               timfcn=timfcn+dum4-dum1
               call timing(dum1,dum2,dum3)
            endif
c
c           --- fout(*,2)=df/d(rhoa), 
c               fout(*,3)=df/d(gaa), 
c               fout(*,4)=df/d(rhob),
c               fout(*,5)=df/d(gbb) , 
c               fout(*,6)=df/d(gab)
c           --- calculate the exchange correlation energy.
            if (calce) then
               exc=exc+sdot(ngtmp,tmpgwt,1,fout,1)
            endif
c           --- now loop over the function pairs and form 
c               the exchange/corr matrix elements
            if(timeit) then
               call timing(dum1,dum2,dum3)
            endif
            call fmkclos(ngtmp,mxgbsiz,minesz,nbf,nnp,becke,lyp,
     $                   tmpgwt,scr,phibar,gradbar,
     $                   queue,tea,fout,phi,grad,dengrada,
     $                   kmat,kmcut)
            if(timeit) then
               call timing(dum4,dum5,dum6)
               timk=timk+dum4-dum1
            endif
c
c
 14         continue 
            ioff=ioff+ng
 15      continue 
c         write(6,*)" Node:",mynodeid," finished with atom ",iatom
         next=nxtval(nproc)+1
 10   continue 
c     

c      if(timeit) then
c         write(4,*) 'node',mynodeid,
c     $        ': time_bf,time_den,time_sqz,time_fcn,time_k',
c     $                  timbf,timden,timsqz,timfcn,timk
c      endif
c
c
      next=nxtval(-nproc)
      return
      end
