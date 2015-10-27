*deck %W%  %G%
CRLM this is a temporary fix so that dxopen looks like dxclos.f
CRLM it may be a good starting place to make open shell gradients work.
CRLM take a look at dxopen_tom for a first attempt by Tom Russo.
      subroutine dxopen(dengrida,dengrada,ga,
     $     fout,itch,phi,grad,hess,scr,slater,becke,lyp,vwn,d,nbf,nnp,
     $     dkay,grdwt,gradwts,ndmat,nat,mxgrd,calce,exc,
     $     dmcut,dencut,calc,nzptrs,tmpgwt,eexch,prteexch,queue,tea,
     $     mxgbsiz,mxgblk,
     $     c,ex,cont,ptprim,noprim,nocont,ptcont,
     $     mxcont,nprim,ntypes,nbtype,ncont,
     $     start,nocart,nobf,maxmom,minmom,mintyp,
     $     nx,ny,nz,xyzgrid,maxl,bigl,charge,grade,accum,ian,ptrad,
     $     radshls,vwts,rnuc,amu,pwtx,rr,radii,akl,rmax,lmax,nomega,
     $     nradial,grdtyp,adjust,minesz)
c***begin prologue     %M%
c***date written       930521   (yymmdd)  
c***revision date      %G%
c   may 13, 1994       rlm at lanl
c      modfying kmtxgts(m511) to do derivatives.
c
c***keywords           k matrix coulomb matrix closed shell
c***author             RUSSO, thomas (lanl)
c***source             %W%   %G%
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
c***end prologue       %M%
      implicit none
c     --- input variables -----
      integer ndmat,nat,mxgrd,nnp,nbf,mxgblk,mxgbsiz
      integer nprim,ntypes,nbtype,ncont,mxcont
      integer rmax,lmax,nomega,nradial
      integer minesz
      real*8 dencut,dmcut
      logical slater,becke,lyp,vwn,calce,prteexch
      character*(*) calc
c     --- input arrays (unmodified) ---
      integer ptprim(nat,ntypes),noprim(nat,ntypes),nocont(nat,ntypes)
      integer ptcont(nat,ntypes),start(nat,ntypes),nocart(0:*)
      integer nobf(ntypes),maxmom(ntypes),minmom(ntypes),mintyp(ntypes)
      integer nx(*),ny(*),nz(*)
      integer ian(nat)
      real*8 c(3,nat),ex(nprim),cont(ncont)
      real*8 xyzgrid(mxgrd,3)
      real*8 d(nnp),exc,eexch,grdwt(mxgrd),gradwts(mxgrd,3,nat)
      logical adjust
      character*(*) grdtyp(nat)
c     --- input arrays (scratch) ---
      integer nzptrs(mxgbsiz)
      integer maxl(nat)
      integer ptrad(rmax,nat),radshls(nat)
      real*8 dengrida(mxgbsiz),dengrada(mxgbsiz,3),
     $       ga(mxgbsiz),
     $       fout(mxgbsiz,*),itch(*),phi(mxgbsiz,nbf),
     $       grad(mxgbsiz,nbf,3),scr(mxgbsiz),hess(mxgbsiz,nbf,6)
      real*8 tmpgwt(mxgbsiz)
      real*8 queue(mxgbsiz,3),tea(mxgbsiz,nbf),accum(nbf,3)
      real*8 vwts(mxgrd),rnuc(nat,nat),amu(nat,nat),pwtx(nat),
     $       rr(nat),radii(nat),akl(nat,nat)
c     --- output arrays ---
      real*8 dkay(nbf,nbf,3)
      real*8 charge(nat)
      real*8 grade(3,nat)
c     --- output variables ---
c     --- scratch arrays ---
c     --- local variables ---
      integer ngtmp,igblk,ioff,ng
      integer rho,zeta,g,h,t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11
      integer x,t,omega,delta,iatom,mu,nu,derivs,i,j,jatom,jtype,ibf
      integer dsq
      integer inp,iout
      integer xyzpow,r,s,tee,rsq
      integer top,ngot,wpadti
      integer bigl,coord
      integer ngb,ngrid,oddblock,phibar,gradbar
      integer hescod(3,3)
      real*8 zero,half,one,two,four
      real*8 sdot
      real*8 slalph,beckb,lypa,lypb,lypc,lypd
      real*8 dum1,dum2,dum3,dum4,dum5,dum6
      real*8 timbf,timden,timsqz,timfcn,timk
      real*8 g1,g2,g3
      logical timeit
      logical fgrad,debug
c
c     ---  used to find a given second derivative in a hessian as 
c          i've laid it out:
c          if you want d2/dx1dx2 you grab hess(.,bf,hescod(x1,x2))
      data hescod/1,4,5,4,2,6,5,6,3/
c
      data timbf,timden,timsqz,timfcn,timk/5*0.0d+00/
      save timbf,timden,timsqz,timfcn,timk
c
      common /io/ inp,iout
c
      parameter (timeit=.false.,debug=.false.)
      parameter (zero=0.0d+00,half=0.5d+00,one=1.0d+00,two=2.0d+00)
      parameter (four=4.0d0)
      parameter (slalph=2.0d0/3.0d0,beckb=0.0042d0,lypa=0.04918d0)
      parameter (lypb=0.132d0,lypc=0.2533d0,lypd=0.349d0)
c
c
      integer mynodeid,nodeid,nproc,nnodes,mdtob,mitob,next,nxtval
      include 'msgtypesf.h'
c
 1010 format(5x,80a1)
c
c      
      mynodeid=nodeid()
      nproc=nnodes()
c
c     --- does the functional involve gradient terms.
      fgrad=.false.
      if(becke.or.lyp) then
         fgrad=.true.
      endif

c     --- core allocation within itch
c         note that the space is used by both bfgrd and the functional
c         routines
c         make sure we can really fit in as much as we think
      call getscm(0,itch,ngot,'dxclos',0)
c
c     --- space for bfgrd.
      s=1
      r=s+mxgbsiz*mxcont
      tee=r+mxgbsiz*mxcont
      rsq=tee+mxgbsiz*mxcont
      xyzpow=rsq+mxgbsiz
      top=xyzpow+3*mxgbsiz*bigl
      if(wpadti(top).gt.ngot) then
         if (mynodeid.eq.0)
     $        write(iout,*) 'need,have',wpadti(top),ngot
         call plnkerr('need more memory in dxclos',201)
      endif
c
c     --- space for gridden. can write over bfgrd space.
      phibar=top
      gradbar=top+nbf
      top=gradbar+nbf
c
c     --- functional scratch space. can write over gridden space.
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
         top=delta+mxgbsiz
      endif
      dsq=top
      top=dsq+nbf*nbf
      if(wpadti(top).gt.ngot) then
         if (mynodeid.eq.0)
     $        write(iout,*) 'need,have',wpadti(top),ngot
         call plnkerr('need more core in dxclos',202)
      endif
c     
c     --- loop over atomic grids. this is the main loop for the
c         quadrature.
      call rzero(charge,nat)
      next=nxtval(nproc)+1
c      write(6,*)" Next is ",next
      do 10 iatom=1,nat
         if (iatom .ne. next) goto 10
c         write(6,*)" processor ",mynodeid," doing ",next
c 
c        generate grid, calc number of grid blocks
         call mkatmg(c,ian,xyzgrid,grdwt,rmax,lmax,nomega,
     $        nradial,
     $        nat,ngrid,mxgrd,vwts,rnuc,amu,pwtx,rr,adjust,radii,
     $        akl,grdtyp(iatom),
     $        radshls,ptrad,.true.,gradwts,iatom)
         ngb=ngrid/mxgbsiz
         oddblock=mod(ngrid,mxgbsiz)
         if (oddblock.ne.0) ngb=ngb+1
         ioff=0
         call rzero(dkay,nbf*nbf*3)
         do 15 igblk=1,ngb
c
c           --- calc size of this block
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
     $           nat,nprim,ntypes,nbtype,nnp,ncont,
     $           start,nbf,nocart,nobf,maxmom,minmom,mintyp,
     $           nx,ny,nz,xyzgrid,mxgrd,
     $           itch(xyzpow),scr,itch(rsq),itch(s),itch(r),
     $           itch(tee),phi,grad,hess,.true.,.true.,mxcont,mxgbsiz,
     $           ng,ioff,maxl)
            if(debug.and. mynodeid.eq.0) then
               write(iout,*)" Grid:"
               do 20 i=1,ng
                  write(iout,*) i,xyzgrid(i,1),xyzgrid(i,2),xyzgrid(i,3)
   20          continue
               do 25 i=1,ng
                  write(iout,*)"  ",i,(hess(i,j,1),j=1,nbf)
                  write(iout,*)"      ",(hess(i,j,2),j=1,nbf)
                  write(iout,*)"      ",(hess(i,j,3),j=1,nbf)
                  write(iout,*)"      ",(hess(i,j,4),j=1,nbf)
                  write(iout,*)"      ",(hess(i,j,5),j=1,nbf)
                  write(iout,*)"      ",(hess(i,j,6),j=1,nbf)
   25          continue
            endif
            if(timeit) then
               call timing(dum4,dum5,dum6)
               timbf=timbf+dum4-dum1
            endif
c
c           --- clear the density and density gradient arrays.
            call rzero(dengrida,mxgbsiz)
            if (fgrad) then
               call rzero(dengrada,mxgbsiz*3)
            endif
c
c           --- form density and gradients on the grid for this atom, 
c               also the gradient invariants.
            if(timeit) then
               call timing(dum1,dum2,dum3)
            endif
            call gridden(nbf,nnp,ng,mxgbsiz,d,phi,grad,
     $                   dengrida,dengrada,minesz,scr,itch(phibar),
     $                   itch(gradbar),dmcut,fgrad)
            if(timeit) then
               call timing(dum4,dum5,dum6)
               timden=timden+dum4-dum1
            endif
c
c           --- now we clean the density, and then use nzident to squeeze out 
c               all those values which are zero.  We also need to squeeze the
c               gradient.  once the density's squeezed, squeeze out the grid
c               weights corresponding to zero density points.
            call vclean(dengrida,dencut,ng)
            ngtmp=ng
            call nzident(dengrida,nzptrs,ngtmp)
            if(timeit) then
               call timing(dum1,dum2,dum3)
            endif
            call gthr(dengrida,dengrida,nzptrs,ngtmp)
            do 35 mu=1,nbf
               call gthr(phi(1,mu),phi(1,mu),nzptrs,ngtmp)
               call gthr(grad(1,mu,1),grad(1,mu,1),nzptrs,ngtmp)
               call gthr(grad(1,mu,2),grad(1,mu,2),nzptrs,ngtmp)
               call gthr(grad(1,mu,3),grad(1,mu,3),nzptrs,ngtmp)
               call gthr(hess(1,mu,1),hess(1,mu,1),nzptrs,ngtmp)
               call gthr(hess(1,mu,2),hess(1,mu,2),nzptrs,ngtmp)
               call gthr(hess(1,mu,3),hess(1,mu,3),nzptrs,ngtmp)
               call gthr(hess(1,mu,4),hess(1,mu,4),nzptrs,ngtmp)
               call gthr(hess(1,mu,5),hess(1,mu,5),nzptrs,ngtmp)
               call gthr(hess(1,mu,6),hess(1,mu,6),nzptrs,ngtmp)
 35         continue 
            if (fgrad) then
               call gthr(dengrada(1,1),dengrada(1,1),nzptrs,ngtmp)
               call gthr(dengrada(1,2),dengrada(1,2),nzptrs,ngtmp)
               call gthr(dengrada(1,3),dengrada(1,3),nzptrs,ngtmp)
               call vmul(ga,dengrada,dengrada,ngtmp)
               call vwxy(ga,ga,dengrada(1,2),dengrada(1,2),+1,
     $                   ngtmp)
               call vwxy(ga,ga,dengrada(1,3),dengrada(1,3),+1,
     $                   ngtmp)
               call vclean(ga,dencut,ngtmp)
            endif
            call gthr(tmpgwt,grdwt(ioff+1),nzptrs,ngtmp)
            if(timeit) then
               call timing(dum4,dum5,dum6) 
               timsqz=timsqz+dum4-dum1
            endif

c
c           --- sum the density to obtain atomic charges ---
            charge(iatom)=charge(iatom)
     $                   +two*sdot(ngtmp,dengrida,1,tmpgwt,1)
c
c           --- form functional's derivatives on the grid.
c               fout(*,1)=f(rhoa,rhob,gaa,gbb,gab)
c               fout(*,2)=df/d(rhoa), fout(*,3)=df/d(gaa), 
c               fout(*,4)=df/d(rhob), fout(*,5)=df/d(gbb)
c               fout(*,6)=df/d(gab)
c               both of the exchange functionals clobber the contents 
c               of fout.
            if(timeit) then
               call timing(dum1,dum2,dum3)
            endif
            derivs=-1
            if (calce) derivs=1
c
c           --- exchange pieces ---
            if (slater) then
               call slaterf(ngtmp,mxgbsiz,dengrida,dengrida,slalph,
     $                      derivs,fout)
            else
               call beckef(ngtmp,mxgbsiz,dengrida,dengrida,ga,ga,
     $                    itch(x),itch(t),derivs,beckb,fout,calc,1)
            endif
            if (calce .and. prteexch) then
               eexch=eexch+sdot(ngtmp,tmpgwt,1,fout,1)
            endif
c
c           --- correlation pieces ---
            if (vwn) then
               call vwnf(ngtmp,mxgbsiz,dengrida,dengrida,calc,derivs,
     $                   fout,itch(rho),itch(x),itch(zeta),itch(g),
     $                   itch(h),itch(t1),itch(t2),itch(t3),itch(t4),
     $                   itch(t5),itch(t6),itch(t7),itch(t8),itch(t9),
     $                   itch(t10),itch(t11))
            else if (lyp) then
               call rzero(fout(1,6),mxgbsiz)
               call lypf(ngtmp,mxgbsiz,dengrida,dengrida,ga,ga,ga,
     $                   itch(omega),itch(delta),itch(t),derivs,lypa,
     $                   lypb,lypc,lypd,fout)
            endif
            if (calce) then
               exc=exc+sdot(ngtmp,tmpgwt,1,fout,1)
            endif
            if(timeit) then
               call timing(dum4,dum5,dum6)
               timfcn=timfcn+dum4-dum1
               call timing(dum1,dum2,dum3)
            endif
c
c           --- do the derivative weight trick
c               the gradwts aren't gathered yet, but we can do that now:
            do 6060 jatom=1,nat
               if (jatom .eq. iatom) goto 6060
               do 6061 coord=1,3
                  call gthr(scr,gradwts(ioff+1,coord,jatom),
     $                      nzptrs,ngtmp)
                  g1=sdot(ngtmp,fout,1,scr,1)
                  grade(coord,jatom)=grade(coord,jatom)+g1
                  grade(coord,iatom)=grade(coord,iatom)-g1
 6061          continue 
 6060       continue 
            if(debug.and. mynodeid.eq.0) then
               write(iout,*)"After deriv weight term, grade:"
               do 6066 jatom=1,nat
                  write(iout,*)jatom,grade(1,jatom),grade(2,jatom),
     $                         grade(3,jatom)
 6066          continue 
            endif
c
c           --- now form the "derivative" of the K-matrix.  Sorta.
            if(timeit) then
               call timing(dum1,dum2,dum3)
            endif
c           --- first term (independent of whether functional is
c               density gradient dependent)
            call vmul(scr,fout(1,2),tmpgwt,ngtmp)
            do 60 nu=1,nbf
               call vmul(tea(1,nu),scr,phi(1,nu),ngtmp)
   60       continue
            do 70 coord=1,3
               call sgemm('t','n',nbf,nbf,ngtmp,one,grad(1,1,coord),
     $                    mxgbsiz,tea,mxgbsiz,one,dkay(1,1,coord),nbf)
   70       continue
c
c           --- term dependent on gradient corrected functionals.
c               recall: X(mu,nu)= phi(nu)Grad(Grad phi(mu))^T
c                               +Grad(phi(mu))Grad(phi(nu))^T
c               and so we need to form:
c               dk(mu,nu)=sum(grid) [w*(2*f3+f6)*X(mu,nu)*grad(rho)
c
c               this gives us two terms for the dk(mu,nu,j)
c
c               term1=sum(grid)(wg*(2*f3+f6)*phi(nu)
c                     *sum(i)[H(grid,hescod(j,i),mu)*dengrad(grid,i)]
c               and
c 
c               term2=sum(grid)(wg*(2*f3+f6)*grad(j)(phi(mu))
c                     *sum(i)[grad(i)phi(nu)*dengrad(grid,i)]
            if (becke .or. lyp) then
               call vmove(scr,fout(1,3),ngtmp)
               call smul(scr,scr,two,ngtmp)
               if (lyp) call vadd(scr,scr,fout(1,6),ngtmp)
               call vmul(scr,scr,tmpgwt,ngtmp)
c
c              first the piece dependent on hessians
               do 110 j=1,3
                  call rzero(tea,mxgbsiz*nbf)
                  do 90 i=1,3
                     do 80 mu=1,nbf
                        call vwxy(tea(1,mu),tea(1,mu),
     $                            hess(1,mu,hescod(j,i)),
     $                            dengrada(1,i),+1,ngtmp)
 80                  continue 
 90               continue 
                  do 100 mu=1,nbf
                     call vmul(tea(1,mu),tea(1,mu),scr,ngtmp)
 100              continue 
                  call sgemm('t','n',nbf,nbf,ngtmp,one,tea,mxgbsiz,phi,
     $                       mxgbsiz,one,dkay(1,1,j),nbf)
 110           continue 
c
c              now the other piece, with only first derivs
               call rzero(tea,mxgbsiz*nbf)
               do 130 i=1,3
                  do 120 nu=1,nbf
                     call vwxy(tea(1,nu),tea(1,nu),grad(1,nu,i),
     $                         dengrada(1,i),+1,ngtmp)
 120              continue 
 130           continue 
               do 140 nu=1,nbf
                  call vmul(tea(1,nu),tea(1,nu),scr,ngtmp)
 140           continue 
               do 150 j=1,3
                  call sgemm('t','n',nbf,nbf,ngtmp,one,grad(1,1,j),
     $                       mxgbsiz,tea,mxgbsiz,one,dkay(1,1,j),nbf)
 150           continue 
            endif
            if(timeit) then
               call timing(dum4,dum5,dum6)
               timk=timk+dum4-dum1
            endif
            ioff=ioff+ng
 15      continue 
         if(debug.and.mynodeid.eq.0) then
 1001       format(5x,'the derivative k-matrix')
 1011       format(5x,'coordinate:',i3)
            write(iout,1001)
            do 1012 coord=1,3
               write(iout,1011) coord
               call matout(dkay(1,1,coord),nbf,nbf,nbf,nbf,iout)
 1012       continue 
         endif
c
c        --- we now have the terms we need to calculate this atom's 
c            contribution to the gradient, 
c            which is -2*sum(mu')sum(nu) P(mu,nu)DK(mu,nu).  
c            the mu' summation is over only he basis functions on the atom 
c            for which the gradient is being taken.
c
c            we do it in two steps: first do the unrestricted nu summation, then
c            do a bunch of restricted mus.  There is a trick: the contribution 
c            to del(a)Exc (derivative w.r.t. center A) from grid A are 
c            calculated by translational invariance, whereas del(b)Exc
c            form grid A are calculated as the sum above.
c
         call trtosq(itch(dsq),d,nbf,nnp)
         call fmaccum(itch(dsq),dkay(1,1,1),accum(1,1),nbf)
         call fmaccum(itch(dsq),dkay(1,1,2),accum(1,2),nbf)
         call fmaccum(itch(dsq),dkay(1,1,3),accum(1,3),nbf)
         do 1066 jatom=1,nat
            if (iatom.eq.jatom) goto 1066
            g1=zero
            g2=zero
            g3=zero
            do 1963 jtype=1,nbtype
               if (noprim(jatom,jtype).gt.0) then
                  do 1776 ibf=1,nocont(jatom,jtype)*nobf(jtype)
                     g1=g1-two*accum(ibf+start(jatom,jtype),1)
                     g2=g2-two*accum(ibf+start(jatom,jtype),2)
                     g3=g3-two*accum(ibf+start(jatom,jtype),3)
 1776             continue 
               endif
 1963       continue 
c
c           --- the two is for double occupancy.
            grade(1,jatom)=grade(1,jatom)+two*g1
            grade(2,jatom)=grade(2,jatom)+two*g2
            grade(3,jatom)=grade(3,jatom)+two*g3
            grade(1,iatom)=grade(1,iatom)-two*g1
            grade(2,iatom)=grade(2,iatom)-two*g2
            grade(3,iatom)=grade(3,iatom)-two*g3
 1066    continue 
         if(debug.and.mynodeid.eq.0) then
            write(iout,*)" After atom ",iatom," contribution, grade:"
            do 1067 jatom=1,nat
               write(iout,*)jatom,grade(1,jatom),grade(2,jatom),
     $                      grade(3,jatom)
 1067       continue 
         endif
c
c
         next=nxtval(nproc)+1
   10 continue 
c      write(6,*) "Doing barrier"
      next=nxtval(-nproc)
c     
c
      if(timeit.and.mynodeid.eq.0) then
         write(iout,*) 'time_bf,time_den,time_sqz,time_fcn,time_k',
     $                  timbf,timden,timsqz,timfcn,timk
      endif
c
c
      return
      end
