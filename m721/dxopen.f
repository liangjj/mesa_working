*deck @(#)dxopen.f	5.2 2/5/95
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
c   may 13, 1994       rlm at lanl
c      modifying kmopen(m511) to do derivatives of exchange-correlation
c      energy.
c***keywords           
c***author             russo, thomas(lanl)
c***source             @(#)dxopen.f	5.2   2/5/95
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
c     --- core allocation within itch
c         note that the same region is used for scratch in both
c         bfgrd and the functional routines.
      call getscm(0,itch,ngot,'kmtxguts',0)
      if(directk) then
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
      endif   
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
c loop over atomic grids
c
      call rzero(charge,ndmat*nat)
      bsetptr=0
      bsgptr=0
      gradptr=0
      denptr=0
      do 10 iatom=1,nat
         ioff=0
         do 15 igblk=1,ngb(iatom)
            ng=gblksiz(igblk,iatom)
            if(directk) then
c              --- calculate the amplitudes for this grid block.
               if(timeit) then
                  call timing(dum1,dum2,dum3)
               endif
               call bfgrd(c,ex,cont,ptprim,noprim,nocont,ptcont,
     $                    nat,nprim,ntypes,nbtype,nnp,ncont,
     $                    start,nbf,nocart,nobf,maxmom,minmom,mintyp,
     $                    nx,ny,nz,xyzgrid(1,1,iatom),
     $                    mxgrd,itch(xyzpow),scr,itch(rsq),itch(s),
     $                    itch(r),itch(tbf),phi,grad,itch(hess),
     $                    dograd,.false.,mxcont,mxgbsiz,ng,ioff,maxl)
               if(timeit) then
                  call timing(dum4,dum5,dum6)
                  timbf=timbf+dum4-dum1
               endif
            else
c
c  read in basis functions on grid
c
               call iosys('read real "basis set on grid"'
     $                  //' from rwf',
     $                  mxgbsiz*nbf,phi,bsetptr,' ')
               bsetptr=bsetptr+mxgbsiz*nbf
               if(dograd) then
                  call iosys('read real "basis set gradients'
     $                     //' on grid" from rwf',
     $                     mxgbsiz*nbf*3,grad,gradptr,' ')
                  gradptr=gradptr+mxgbsiz*nbf*3
               endif
            endif
c
c We save "open" and "closed" densities, not "alpha" and "beta"
c this is because of the add after the density formation
c
            if (.not. dorhoup) then
               call rzero(dengrida,mxgbsiz)
               call rzero(dengridb,mxgbsiz)
               if (becke .or. lyp) then
                  call rzero(dengrada,mxgbsiz*3)
                  call rzero(dengradb,mxgbsiz*3)
               endif
            else
               call iosys('read real "old densities" from rwf',
     $              mxgbsiz,dengrida,denptr,' ')
               call iosys('read real "old densities" from rwf',
     $              mxgbsiz,dengridb,denptr+mxgbsiz,' ')
               if (becke .or. lyp) then
                  call iosys('read real "old gradients" from rwf',
     $                 mxgbsiz*3,dengrada,denptr*3,' ')
                  call iosys('read real "old gradients" from rwf',
     $                 mxgbsiz*3,dengradb,denptr*3+3*mxgbsiz,' ')
               endif
            endif
c
c form density and gradients on the grid for this atom, also the gradient
c invariants.
c
            kl=0
            if(timeit) then
               call timing(dum1,dum2,dum3)
            endif
            do 20 mu=1,nbf
               do 30 nu=1,mu
                  kl=kl+1
                  if (abs(d(kl,1)).lt.dmcut .and.
     $                 abs(d(kl,2)).lt.dmcut) goto 30
c     
c     dengrid(i)=sum over mu,nu P(mu,nu)*phi(i,mu)*phi(i,nu)
c     dengrida=alpha, dengridb=beta.  Alpha is beta+open, so first
c     form open/closed.  Put closed in dengridb, open in dengrida, fix
c     later.
                  if (mu.eq.nu) then
                     call vmul(scr,phi(1,mu),phi(1,nu),ng)
                     call vwxs(dengridb,dengridb,scr,d(kl,1),+1,
     $                    ng)
                     call vwxs(dengrida,dengrida,scr,d(kl,2),+1,
     $                    ng)
                     if (becke .or. lyp) then
                        foo1=two*d(kl,1)
                        foo2=two*d(kl,2)
                        do 31 coord=1,3
                           call vmul(scr,phi(1,mu),grad(1,nu,coord),
     $                          ng)
                           call vwxs(dengradb(1,coord),
     $                          dengradb(1,coord),scr,foo1,+1,ng)
                           call vwxs(dengrada(1,coord),
     $                          dengrada(1,coord),scr,foo2,+1,ng)
 31                     continue 
                     endif
                  else
                     foo1=two*d(kl,1)
                     foo2=two*d(kl,2)
                     call vmul(scr,phi(1,mu),phi(1,nu),ng)
                     call vwxs(dengridb,dengridb,scr,foo1,+1,ng)
                     call vwxs(dengrida,dengrida,scr,foo2,+1,ng)
                     if (becke .or. lyp) then
                        do 32 coord=1,3
                           call vmul(scr,phi(1,mu),grad(1,nu,coord),
     $                          ng)
                           call vwxy(scr,scr,phi(1,nu),grad(1,mu,coord),
     $                          +1,ng)
                           call vwxs(dengradb(1,coord),
     $                          dengradb(1,coord),scr,foo1,+1,ng)
                           call vwxs(dengrada(1,coord),
     $                          dengrada(1,coord),scr,foo2,+1,ng)
 32                     continue 
                     endif
                  endif
 30            continue 
 20         continue 
            if(timeit) then
               call timing(dum4,dum5,dum6)
               timden=timden+dum4-dum1
            endif
            call iosys('write real "old densities" on rwf',
     $           mxgbsiz,dengrida,denptr,' ')
            call iosys('write real "old densities" on rwf',
     $           mxgbsiz,dengridb,denptr+mxgbsiz,' ')
            if (becke .or. lyp) then
               call iosys('write real "old gradients" on rwf',
     $              mxgbsiz*3,dengrada,denptr*3,' ')
               call iosys('write real "old gradients" on rwf',
     $              mxgbsiz*3,dengradb,denptr*3+mxgbsiz*3,' ')
            endif
            denptr=denptr+mxgbsiz*2
c
c now combine open with closed to get alpha
c
            call vadd(dengrida,dengrida,dengridb,ng)
            if (becke .or. lyp) then
               call vadd(dengrada(1,1),dengrada(1,1),dengradb(1,1),
     $              ng)
               call vadd(dengrada(1,2),dengrada(1,2),dengradb(1,2),
     $              ng)
               call vadd(dengrada(1,3),dengrada(1,3),dengradb(1,3),
     $              ng)
            endif
c
c clean based on total density
c
            call vclean(dengrida,dencut,ng)
            call vclean(dengridb,dencut,ng)
c form total (alpha+beta) density using tmpgwt as a temp store area
            call vadd(tmpgwt,dengrida,dengridb,ng)
            ngtmp=ng 
            call nzident(tmpgwt,nzptrs,ngtmp)
            if(timeit) then
               call timing(dum1,dum2,dum3)
            endif
            call gthr(dengrida,dengrida,nzptrs,ngtmp)
            call gthr(dengridb,dengridb,nzptrs,ngtmp)
            call gthr(tmpgwt,grdwt(1+ioff,iatom),nzptrs,ngtmp)
            if (becke.or.lyp)then
               call gthr(dengrada,dengrada,nzptrs,ngtmp)
               call gthr(dengrada(1,2),dengrada(1,2),nzptrs,ngtmp)
               call gthr(dengrada(1,3),dengrada(1,3),nzptrs,ngtmp)
               call gthr(dengradb,dengradb,nzptrs,ngtmp)
               call gthr(dengradb(1,2),dengradb(1,2),nzptrs,ngtmp)
               call gthr(dengradb(1,3),dengradb(1,3),nzptrs,ngtmp)
            endif
            do 35 mu=1,nbf
               call gthr(phi(1,mu),phi(1,mu),nzptrs,ngtmp)
               if (becke.or.lyp)then
                  call gthr(grad(1,mu,1),grad(1,mu,1),nzptrs,ngtmp)
                  call gthr(grad(1,mu,2),grad(1,mu,2),nzptrs,ngtmp)
                  call gthr(grad(1,mu,3),grad(1,mu,3),nzptrs,ngtmp)
               endif
 35         continue 
            if (becke .or. lyp) then
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
c at this point we have the alpha and beta densities and their gradients.
c Now we need to construct the k matrices.
c
            if(timeit) then
               call timing(dum1,dum2,dum3)
            endif
            derivs=-1
            if (calce) derivs=1
            if (lyp) call rzero(fout(1,6),mxgbsiz)
            if (slater) then
               call slaterf(ngtmp,mxgbsiz,dengrida,dengridb,slalph,
     $              derivs,fout)
            else
c
c use 'uhf' because we really are doing a general beckef call...may be time
c to remove that silliness?
c     
               call beckef(ngtmp,mxgbsiz,dengrida,dengridb,ga,gb,
     $              itch(x),itch(t),derivs,beckb,fout,'uhf',1)
            endif
            if (calce .and. prteexch) then
               eexch=eexch+sdot(ngtmp,fout(1,1),1,tmpgwt,1)
            endif
            if (vwn) then
               call vwnf(ngtmp,mxgbsiz,dengrida,dengridb,calc,derivs,
     $              fout,itch(rho),itch(x),itch(zeta),itch(g),
     $              itch(h),itch(t1),itch(t2),itch(t3),itch(t4),
     $              itch(t5),itch(t6),itch(t7),itch(t8),itch(t9),
     $              itch(t10),itch(t11))
            else if (lyp) then
               call lypf(ngtmp,mxgbsiz,dengrida,dengridb,ga,gb,gab,
     $              itch(omega),itch(delta),itch(t),derivs,lypa,lypb,
     $              lypc,lypd,fout)
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
c part common to all functionals
c
            call vmul(scr,fout(1,2),tmpgwt,ngtmp)
            do 60 mu=1,nbf
               call vmul(tea(1,mu),scr,phi(1,mu),ngtmp)
 60         continue 
            call sgemm('t','n',nbf,nbf,ngtmp,one,phi,mxgbsiz,tea,
     $           mxgbsiz,zero,kay,nbf)
            call sqtotr(ktmp,kay,nbf,nnp)
            call vadd(kmat(1,1),kmat(1,1),ktmp,nnp)
            call vmul(scr,fout(1,4),tmpgwt,ngtmp)
            do 70 mu=1,nbf
               call vmul(tea(1,mu),scr,phi(1,mu),ngtmp)
 70         continue 
            call sgemm('t','n',nbf,nbf,ngtmp,one,phi,mxgbsiz,tea,
     $           mxgbsiz,zero,kay,nbf)
            call sqtotr(ktmp,kay,nbf,nnp)
            call vadd(kmat(1,2),kmat(1,2),ktmp,nnp)
c
c gradient-dependent bit
c
            if (becke .or. lyp) then
c     
c alpha part
c
               call vmove(scr,fout(1,3),ngtmp)
               call smul(scr,scr,two,ngtmp)
               call vmul(scr,scr,tmpgwt,ngtmp)
               if  (lyp) then 
                  call  vmove(scr2,fout(1,6),ngtmp)
                  call vmul(scr2,scr2,tmpgwt,ngtmp)
               endif             
               do 80 i=1,3
                  call vmul(queue(1,i),scr,dengrada(1,i),ngtmp)
                  if (lyp)
     $                 call vwxy(queue(1,i),queue(1,i),scr2,
     $                 dengradb(1,i),+1,ngtmp)
 80            continue 
               do 90 mu=1,nbf
                  call vmul(tea(1,mu),queue(1,1),grad(1,1,mu),ngtmp)
                  call vwxy(tea(1,mu),tea(1,mu),queue(1,2),grad(1,mu,2),
     $                 +1,ngtmp)
                  call vwxy(tea(1,mu),tea(1,mu),queue(1,3),grad(1,mu,3),
     $                 +1,ngtmp)
 90            continue 
               call sgemm('t','n',nbf,nbf,ngtmp,one,phi,mxgbsiz,tea,
     $              mxgbsiz,zero,kay,nbf)
               do 100 mu=1,nbf
                  kay(mu,mu)=2*kay(mu,mu)
 100           continue 
               do 110 mu=1,nbf
                  do 120 nu=1,mu-1
                     kay(mu,nu)=kay(mu,nu)+kay(nu,mu)
                     kay(nu,mu)=kay(mu,nu)
 120              continue 
 110           continue 
               call sqtotr(ktmp,kay,nbf,nnp)
               call vadd(kmat(1,1),kmat(1,1),ktmp,nnp)
c     
c     beta part
c
               call vmove(scr,fout(1,5),ngtmp)
               call smul(scr,scr,two,ngtmp)
               call vmul(scr,scr,tmpgwt,ngtmp)
               if  (lyp) then
                  call vmove(scr2,fout(1,6),ngtmp)
                  call vmul(scr2,scr2,tmpgwt,ngtmp)
               endif
               do 180 i=1,3
                  call vmul(queue(1,i),scr,dengradb(1,i),ngtmp)
                  if (lyp)
     $                 call vwxy(queue(1,i),queue(1,i),scr2,
     $                 dengrada(1,i),+1,ngtmp)
 180           continue 
               do 190 mu=1,nbf
                  call vmul(tea(1,mu),queue(1,1),grad(1,1,mu),ngtmp)
                  call vwxy(tea(1,mu),tea(1,mu),queue(1,2),grad(1,2,mu),
     $                 +1,ngtmp)
                  call vwxy(tea(1,mu),tea(1,mu),queue(1,3),grad(1,3,mu),
     $                 +1,ngtmp)
 190           continue 
               call sgemm('t','n',nbf,nbf,ngtmp,one,phi,mxgbsiz,tea,
     $              mxgbsiz,zero,kay,nbf)
               do 1100 mu=1,nbf
                  kay(mu,mu)=2*kay(mu,mu)
 1100          continue 
               do 1110 mu=1,nbf
                  do 1120 nu=1,mu-1
                     kay(mu,nu)=kay(mu,nu)+kay(nu,mu)
                     kay(nu,mu)=kay(mu,nu)
 1120             continue 
 1110          continue 
               call sqtotr(ktmp,kay,nbf,nnp)
               call vadd(kmat(1,2),kmat(1,2),ktmp,nnp)
            endif
            if(timeit) then
               call timing(dum4,dum5,dum6)
               timk=timk+dum4-dum1
            endif
            ioff=ioff+ng
 15      continue 
 10   continue 
c
c
      if(timeit) then
         write(iout,*) 'time_bf,time_den,time_sqz,time_fcn,time_k',
     $                  timbf,timden,timsqz,timfcn,timk
      endif
c
c
      return
      end
