*deck @(#)kmopen.f	5.1  11/28/95
      subroutine kmopen(dengrida,dengridb,dengrada,dengradb,ga,gb,gab,
     $     fout,itch,phi,grad,scr,scr2,slater,becke,lyp,vwn,d,
     $     nbf,nnp,kmat,nexch,grdwt,ndmat,nat,mxgrd,calce,exc,
     $     dmcut,dencut,calc,nzptrs,tmpgwt,eexch,
     $     prteexch,queue,tea,dorhoup,mxgbsiz,
     $     c,ex,cont,ptprim,noprim,nocont,ptcont,
     $     mxcont,nprim,ntypes,nbtype,ncont,
     $     start,nocart,nobf,maxmom,minmom,mintyp,
     $     nx,ny,nz,ian,xyzgrid,charge,maxl,dograd,bigl,
     $     vwts,rnuc,amu,pwtx,rr,radii,akl,ptrad,
     $     rmax,lmax,nomega,nradial,grdtyp,adjust,minesz,phibar,
     $     gradbar,shlmax,coeff,bfcut,kmcut,grdfil)
c***begin prologue     kmopen.f
c***date written       yymmdd  
c***revision date      11/28/95      
c
c***keywords           
c***author             
c***source             @(#)kmopen.f	5.1   11/28/95
c***purpose            
c***description
c     
c    
c
c***references
c
c***routines called
c
c***end prologue       kmopen.f
      implicit none
c     --- input variables -----
      integer nat,mxgrd,ndmat,nexch,nnp,nbf,mxgbsiz
      integer nprim,ntypes,nbtype,ncont,mxcont,bigl
      integer rmax,lmax,nomega,nradial
      integer ngrid,minesz
      real*8 dencut,dmcut,bfcut,kmcut
      logical slater,becke,lyp,vwn,calce,prteexch,dorhoup
      logical dograd,adjust
      character*(*) calc
      character*4 grdfil
c     --- input arrays (unmodified) ---
      integer ian(nat)
      integer ptprim(nat,ntypes),noprim(nat,ntypes),nocont(nat,ntypes)
      integer ptcont(nat,ntypes),start(nat,ntypes),nocart(0:*)
      integer nobf(ntypes),maxmom(ntypes),minmom(ntypes),mintyp(ntypes)
      integer nx(*),ny(*),nz(*)
      integer shlmax(2)
      real*8 c(3,nat),ex(nprim),cont(ncont),coeff(nbf,nbf)
      real*8 d(nnp,ndmat)
      character*(*) grdtyp(nat)
c     --- input arrays (scratch) ---
      integer nzptrs(mxgbsiz)
      integer maxl(nat)
      integer ptrad(rmax)
      real*8 dengrida(mxgbsiz),dengridb(mxgbsiz),dengrada(mxgbsiz,3),
     $       dengradb(mxgbsiz,3),ga(mxgbsiz),gb(mxgbsiz),gab(mxgbsiz),
     $       fout(mxgbsiz,*),itch(*),phi(mxgbsiz,nbf),
     $       grad(mxgbsiz,nbf,3),scr(mxgbsiz),scr2(mxgbsiz)
      real*8 phibar(nbf),gradbar(nbf)
      real*8 grdwt(mxgrd)
      real*8 tmpgwt(mxgbsiz)
      real*8 queue(mxgbsiz,3),tea(mxgbsiz,nbf)
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
      integer rho,zeta,g,h,t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11
      integer inp,iout
      integer top,ngot,wpadti
      integer ngtmp,denptr,ioff,ng,igblk
      integer xyzpow,r,s,rsq,tbf,hess
      integer radshls,oddblock,ngb
      logical dogrdio
      real*8 slalph,beckb,lypa,lypb,lypc,lypd
      real*8 zero,half,one,two,four
      real*8 sdot
      real*8 dum1,dum2,dum3,dum4,dum5,dum6
      real*8 timbf,timden,timsqz,timfcn,timk
      logical timeit
c
      data timbf,timden,timsqz,timfcn,timk/5*0.0d+00/
      save timbf,timden,timsqz,timfcn,timk
c
      common /io/ inp,iout
c
      parameter (timeit=.true.)
      parameter (dogrdio=.false.)
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
 1010 format(1x,80a1)
c
c
c 
      mynodeid=nodeid()
      nproc=nnodes()
c
c     --- core allocation within itch
c         note that the same region is used for scratch in both
c         bfgrd and the functional routines.
      call getscm(0,itch,ngot,'kmopen',0)
      s=1
      r=s+mxgbsiz*mxcont
      tbf=r
      hess=tbf
      rsq=r+mxgbsiz*mxcont
      xyzpow=rsq+mxgbsiz
      top=xyzpow+3*mxgbsiz*bigl
      if(wpadti(top).gt.ngot) then
         write(iout,*) 'need,have',wpadti(top),ngot
         call lnkerr('need more memory in kmopen')
      endif
c
c     --- space for functionals.
      t=0
      if (slater) then
c        needs nothing
         top=1
      else if (becke) then
         x=1
         t=x+mxgbsiz
         top=t+4*mxgbsiz
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
         write(iout,*) 'need,have',wpadti(top),ngot
         call lnkerr('need more core in kmopen')
      endif
c
c     --- loop over atomic grids and form the exchange operator.
      if(dogrdio) then
c        can use scratch area ptrad to store grid sizes
         call iosys('read integer "external atomic grid size" from rwf',
     $              nat,ptrad,0,' ') 
      endif
      call rzero(charge,ndmat*nat)
      denptr=0
      next=nxtval(nproc)+1
      do 10 iatom=1,nat
         if (iatom .ne. next) goto 10
c
c        --- generate grid, calc number of grid blocks
         if(dogrdio) then
            ngrid=ptrad(iatom)
            call iosys('read real "external grid" from '//grdfil(1:3),
     $                 mxgrd*3,xyzgrid,mxgrd*3*(iatom-1),' ')
            call iosys('read real "external grid weights" from '
     $                 //grdfil(1:3),mxgrd,grdwt,mxgrd*(iatom-1),' ')
         else
            call mkatmg(c,ian,xyzgrid,grdwt,rmax,lmax,nomega,
     $                  nradial,
     $                  nat,ngrid,mxgrd,vwts,rnuc,amu,pwtx,rr,adjust,
     $                  radii,akl,grdtyp(iatom),
     $                  radshls,ptrad,.false.,vwts,iatom)
         endif
c        --- loop over the large grid blocks.
         ngb=ngrid/mxgbsiz
         oddblock=mod(ngrid,mxgbsiz)
         if (oddblock.ne.0) ngb=ngb+1
         ioff=0
         do 15 igblk=1,ngb
c
c           --- calc size of this block
            if (igblk .ne. ngb .or. oddblock.eq.0 ) then
               ng=mxgbsiz
            else
               ng=oddblock
            endif
c           --- calculate the amplitudes for this grid block.
            if(timeit) then
               call timing(dum1,dum2,dum3)
            endif
            call bfgrd(c,ex,cont,ptprim,noprim,nocont,ptcont,
     $                 nat,nprim,ntypes,nbtype,nnp,ncont,
     $                 start,nbf,nocart,nobf,maxmom,minmom,mintyp,
     $                 nx,ny,nz,xyzgrid,
     $                 mxgrd,itch(xyzpow),scr,itch(rsq),itch(s),
     $                 itch(r),itch(tbf),phi,grad,itch(hess),
     $                 dograd,.false.,mxcont,mxgbsiz,ng,ioff,maxl,
     $                 bfcut)
            if(timeit) then
               call timing(dum4,dum5,dum6)
               timbf=timbf+dum4-dum1
            endif
c
c           ---  we save "open" and "closed" densities, not "alpha" and "beta"
c                this is because of the add after the density formation
            if (.not. dorhoup) then
               call rzero(dengrida,mxgbsiz)
               call rzero(dengridb,mxgbsiz)
               if (becke .or. lyp) then
                  call rzero(dengrada,mxgbsiz*3)
                  call rzero(dengradb,mxgbsiz*3)
               endif
            else
               call plnkerr('Yo! who changed dorhoup?',1103)
               call iosys('read real "open density" from rwf',
     $                     mxgbsiz,dengrida,denptr,' ')
               call iosys('read real "open density" from rwf',
     $                     mxgbsiz,dengridb,denptr+mxgbsiz,' ')
               if (becke .or. lyp) then
                  call iosys('read real "open density gradient"'
     $                    //' from rwf',mxgbsiz*3,dengrada,denptr*3,' ')
                  call iosys('read real "open density gradient"'
     $                     //' from rwf',mxgbsiz*3,dengradb,
     $                                   denptr*3+3*mxgbsiz,' ')
               endif
            endif
c
c           --- form density and gradients on the grid for this atom, 
c               also the gradient invariants.
c               dengrida=alpha, dengridb=beta. 
c                  alpha is beta+open, so first from open/closed.
c                  put closed in dengridb, open in dengrida, 
c                  fix later.
            if(timeit) then
               call timing(dum1,dum2,dum3)
            endif
            call moden(nbf,shlmax(1),ng,mxgbsiz,coeff,phi,grad,
     $                 dengridb,dengradb,minesz,scr,scr2,
     $                 phibar,gradbar,bfcut,dograd)
            call moden(nbf,shlmax(2)-shlmax(1),ng,mxgbsiz,
     $                 coeff(1,shlmax(1)+1),phi,grad,
     $                 dengrida,dengrada,minesz,scr,scr2,
     $                 phibar,gradbar,bfcut,dograd)
c            call gridden(nbf,nnp,ng,mxgbsiz,d(1,1),phi,grad,
c     $                   dengridb,dengradb,minesz,scr,phibar,
c     $                   gradbar,dmcut,dograd)
c            call gridden(nbf,nnp,ng,mxgbsiz,d(1,2),phi,grad,
c     $                   dengrida,dengrada,minesz,scr,phibar,
c     $                   gradbar,dmcut,dograd)
            if(timeit) then
               call timing(dum4,dum5,dum6)
               timden=timden+dum4-dum1
            endif
c$$$            call iosys('write real "open density" on rwf',
c$$$     $                  mxgbsiz,dengrida,denptr,' ')
c$$$            call iosys('write real "open density" on rwf',
c$$$     $                  mxgbsiz,dengridb,denptr+mxgbsiz,' ')
c$$$            if (becke .or. lyp) then
c$$$               call iosys('write real "open density gradient" on rwf',
c$$$     $                     mxgbsiz*3,dengrada,denptr*3,' ')
c$$$               call iosys('write real "open density gradient" on rwf',
c$$$     $                     mxgbsiz*3,dengradb,denptr*3+mxgbsiz*3,' ')
c$$$            endif
            denptr=denptr+mxgbsiz*2
c
c           --- now combine open with closed to get alpha
            call vadd(dengrida,dengrida,dengridb,ng)
            if (becke .or. lyp) then
               call vadd(dengrada(1,1),dengrada(1,1),dengradb(1,1),
     $                   ng)
               call vadd(dengrada(1,2),dengrada(1,2),dengradb(1,2),
     $                   ng)
               call vadd(dengrada(1,3),dengrada(1,3),dengradb(1,3),
     $                   ng)
            endif
c
c           ---  clean based on total density
            if(timeit) then
               call timing(dum1,dum2,dum3)
            endif
            call sqzopn(ng,ngtmp,mxgbsiz,nbf,phi,grad,dengrida,
     $                  dengrada,dengridb,dengradb,grdwt(1+ioff),
     $                  tmpgwt,nzptrs,dencut,becke,lyp) 
c
c           --- form the density invariants.
            if (becke .or. lyp) then
               call vmul(gb,dengradb,dengradb,ngtmp)
               call vwxy(gb,gb,dengradb(1,2),dengradb(1,2),+1,
     $                   ngtmp)
               call vwxy(gb,gb,dengradb(1,3),dengradb(1,3),+1,
     $                   ngtmp)
               call vmul(ga,dengrada,dengrada,ngtmp)
               call vwxy(ga,ga,dengrada(1,2),dengrada(1,2),+1,
     $                   ngtmp)
               call vwxy(ga,ga,dengrada(1,3),dengrada(1,3),+1,
     $                   ngtmp)
               if (lyp) then
                  call vmul(gab,dengrada,dengradb,ngtmp)
                  call vwxy(gab,gab,dengrada(1,2),dengradb(1,2),+1,
     $                      ngtmp)
                  call vwxy(gab,gab,dengrada(1,3),dengradb(1,3),+1,
     $                      ngtmp)
                  call vclean(gab,dencut,ngtmp)
               endif
               call vclean(ga,dencut,ngtmp)
               call vclean(gb,dencut,ngtmp)
            endif
            if(timeit) then
               call timing(dum4,dum5,dum6)
               timsqz=timsqz+dum4-dum1
            endif
c
c           --- sum to get atomic charge ---
            charge(iatom,1)=charge(iatom,1)
     $                     +sdot(ngtmp,tmpgwt,1,dengrida,1)
            charge(iatom,2)=charge(iatom,2)
     $                     +sdot(ngtmp,tmpgwt,1,dengridb,1)
c
c           --- at this point we have the alpha and beta densities and their 
c               gradients. now we construct the k matrices.
            if(timeit) then
               call timing(dum1,dum2,dum3)
            endif
            derivs=-1
            if (calce) derivs=1
            if (lyp) call rzero(fout(1,6),mxgbsiz)
            if (slater) then
               call slaterf(ngtmp,mxgbsiz,dengrida,dengridb,slalph,
     $                     derivs,fout)
            else
c
c              --- use 'uhf' because we really are doing a general beckef 
c                  call... may be time to remove that silliness?
               call beckef(ngtmp,mxgbsiz,dengrida,dengridb,ga,gb,
     $                    itch(x),itch(t),derivs,beckb,fout,'uhf',1)
            endif
            if (calce .and. prteexch) then
               eexch=eexch+sdot(ngtmp,fout(1,1),1,tmpgwt,1)
            endif
            if (vwn) then
               call vwnf(ngtmp,mxgbsiz,dengrida,dengridb,calc,derivs,
     $                   fout,itch(rho),itch(x),itch(zeta),itch(g),
     $                   itch(h),itch(t1),itch(t2),itch(t3),itch(t4),
     $                   itch(t5),itch(t6),itch(t7),itch(t8),itch(t9),
     $                   itch(t10),itch(t11))
            else if (lyp) then
               call lypf(ngtmp,mxgbsiz,dengrida,dengridb,ga,gb,gab,
     $                   itch(omega),itch(delta),itch(t),derivs,lypa,
     $                   lypb,lypc,lypd,fout)
            endif
            if (calce) then
               exc=exc+sdot(ngtmp,fout(1,1),1,tmpgwt,1)
            endif
            if(timeit) then
               call timing(dum4,dum5,dum6)
               timfcn=timfcn+dum4-dum1
               call timing(dum1,dum2,dum3)
            endif
c
c           --- form the exchange matrix.
            call fmkopen(ngtmp,mxgbsiz,minesz,nbf,nnp,becke,lyp,
     $                   tmpgwt,scr,scr2,
     $                   queue,tea,fout,phi,grad,phibar,gradbar,
     $                   dengrada,dengradb,kmat,kmcut)
            if(timeit) then
               call timing(dum4,dum5,dum6)
               timk=timk+dum4-dum1
            endif
c
c
            ioff=ioff+ng
 15      continue 
         next=nxtval(nproc)+1
 10   continue 
c
c
c      if(timeit) then
c         write(iout,*) 'time_bf,time_den,time_sqz,time_fcn,time_k',
c     $                  timbf,timden,timsqz,timfcn,timk
c      endif
c
c

c BARRIER!
      next=nxtval(-nproc)

      return
      end
