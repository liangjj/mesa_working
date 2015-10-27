*deck %W% %G%
      subroutine gofish(itch,d,nbf,nnp,jmat,ncoul,
     $     ndmat,nat,mxgrd,pmxgrd,
     $     dmcut,dencut,bfcut,kmcut,
     $     mxgbsiz,aitch,ian,
     $     c,ex,cont,ptprim,noprim,nocont,ptcont,
     $     mxcont,nprim,ntypes,nbtype,ncont,
     $     start,nocart,nobf,maxmom,minmom,mintyp,
     $     nx,ny,nz,xyzgrid,grdwts,charge,maxl,bigl,ops,
     $     left,vwts,rnuc,amu,pwtx,rr,radii,akl,ptrad,rmax,lmax,
     $     nomega,nradial,grdtyp,adjust,minesz,vlmax,vradial,
     $     vncrule,nlm,rpts,vlm,y2,zdiag,zoffdiag,rmdiag,grdfil)
c***begin prologue     %M%
c***date written       940304      (yymmdd)  
c***revision date      %G%
c
c***keywords           poisson, coulomb, potential, density
c***author             martin, richard(lanl) 
c***source             %W% %G%
c***purpose            generates the coulomb potential from the density
c***description        
c                      solves poisson equation for v, given rho
c                         (del**2) v = rho
c     
c                      this is accomplished by projecting the total
c                      density into a set of atomic pieces using the
c                      functional decomposition of Becke, thereby reducing
c                      the problem to a series of atomic poisson problems.  
c
c                      for each atom, the density is decomposed into
c                      spherical harmonic components and a radial equation
c                      is solved for the potential originating from that
c                      component. the radial equation is converted into an
c                      integral equation using the appropriate green's 
c                      function, and solved via newton-cotes quadrature.
c                      
c                      it should be noted that the grid used to solve the
c                      atomic problem may be  different from that used 
c                      to represent the resulting potential.
c
c***references
c
c***routines called
c
c***end prologue       %M%
      implicit none
c     --- input variables ---
      integer nbf,nnp,ncoul,ndmat,nat,mxgrd,pmxgrd,mxgbsiz
      integer nprim,ntypes,nbtype,ncont,mxcont
      integer rmax,lmax,nradial,bigl
      integer left,minesz,nlm
      integer vlmax,vradial,vncrule
      logical zdiag,zoffdiag,rmdiag
      real*8 dencut,dmcut,bfcut,kmcut
c     --- input arrays (unmodified) ---
      character*(*) ops
      character*(*) grdtyp(nat)
      character*4 grdfil
      integer ian(nat)
      integer ptprim(nat,ntypes),noprim(nat,ntypes),nocont(nat,ntypes)
      integer ptcont(nat,ntypes),start(nat,ntypes),nocart(0:*)
      integer nobf(ntypes),maxmom(ntypes),minmom(ntypes),mintyp(ntypes)
      integer nx(*),ny(*),nz(*)
      real*8 d(nnp,ndmat)
      real*8 c(3,nat),ex(nprim),cont(ncont)
      real*8 xyzgrid(mxgrd,3),grdwts(mxgrd)
c     --- input arrays (scratch) ---
      real*8 itch(left),vwts(mxgrd)
      real*8 rnuc(nat,nat),amu(nat,nat),pwtx(nat),rr(nat)
      real*8 akl(nat,nat),ptrad(*),rpts(vradial,nat)
      real*8 vlm(vradial,nlm,nat),y2(vradial,nlm,nat)
      integer aitch(*)
      integer maxl(nat)
c     --- output arrays ---
      real*8 jmat(nnp,ncoul)
      real*8 charge(nat,ndmat)
c     --- output variables ---
c     --- scratch arrays ---
      logical didvor(1000)
c     --- local variables ---
      logical adjust,rtox,dogrdio,method1
      integer vnomega,jngrid
      integer wpadti,iadtwp
      integer nomega,nr,ngrid,lebord
      integer v,radii,angpts,angwt,top,top0,rwt,grid,gridwts
      integer phibf,scrgrd,rsq,r,s,t,grad,hess,xyzpow,tea,vptrad
      integer dengrid,dengrad
      integer tvwts,twts
      integer ylm,ptlm,ctheta,phi,cphi,sphi,plm,scr,tmp
      integer rholm
      integer ind,fr
      integer i,ioff,joff,iatom,jatom,ngj,ig
      integer lm,l,m
      integer inp,iout
      integer k,jacob
      integer nrings
      integer npring
      integer mxrings,mxcotes
      integer nwtot
      integer j,y
      integer d2
      integer int0f,int0b
      integer phibar,gradbar
      integer ngb,jgblk,radshls
      integer fout,queue
      integer grdsiz
      integer itype,jtype,jtypmx,nconti,ncontj,lenblk
      integer oddblock
      logical debug,dumpfit
      logical timeit,called
      real*8 ylmerr,zero,one
      real*8 test(3)
      real*8 dum1,dum2,dum3,dum4,dum5,dum6
      real*8 timpre,timagrd,timbf,timden,tim2ylm,timrsolv
      real*8 timatg,tim2grd,timfmk
      real*8 sdot
c
      parameter (debug=.false.,timeit=.false.,dumpfit=.false.)
      parameter (rtox=.true.,dogrdio=.true.)
      parameter (mxrings=100,mxcotes=9)
      parameter (zero=0.0d+00,one=1.0d+00)
      parameter (method1=.false.)
      data timpre,timagrd,timbf,timden,tim2ylm,timrsolv/6*0.0d0/
      data timatg,tim2grd,timfmk/3*0.0d0/
      data called/.false./
      data didvor/1000*.false./
      save called

C TCGMSG
c      
      integer nodeid,mynodeid,mdtob,nxtask,nnodes,nproc,next
      integer nxtval
      integer ierr,stderr
      include 'msgtypesf.h'
      logical ispar
      common /tcgmesa/ ispar
c
      common/io/inp,iout
c      
      integer angsiz
      real*8 rhomax,rhomx 
c
      mynodeid=nodeid()
      nproc=nnodes()
      ierr=stderr()
      timden=0.0
      tim2grd=0.0
      timfmk=0.0
c
c     --- allocate core into itch
      call rzero(charge,nat*ndmat)
c
c     --- generate the angular grid for decomposition.
c            lebedev angular points and weights on the 
c            unit sphere.
      lebord=2*vlmax+1
      vnomega=angsiz(lebord)
c     --- allocate some core. these first few will contain the radial
c         points, the potential, and the spline coefficients.
      angpts=1
      angwt=angpts+3*vnomega
      top0=angwt+3*vnomega
      if (top0 .gt. left) then
         write(iout,*)'not enough core for lebedev !'
         write(iout,*)'have ',left,' need ',top0
         call lnkerr('m511: gofish')
      endif
      if(timeit) then
         call timing(dum1,dum2,dum3)
      endif
      call sphere(vnomega,itch(angpts),itch(angwt),lebord,vnomega)
      if(debug) then
         write(iout,*) 'angtps',(itch(angpts+k),k=0,3*vnomega-1)
         write(iout,*) 'angwt',(itch(angwt+k),k=0,vnomega-1)
      endif
c
c     --- generate ylm on the lebedev grid
      ylm=top0
      ptlm=wpadti(ylm+vnomega*nlm)
      ctheta=iadtwp(ptlm+(vlmax+1)*(2*vlmax+1))
      phi=ctheta+vnomega
      cphi=phi+vnomega
      sphi=cphi+vnomega
      plm=sphi+vnomega
      scr=plm+vnomega*(vlmax+1)
      top0=scr+vnomega
      if (top0 .gt. left) then
         write(iout,*)'not enough core for vylm !'
         write(iout,*)'have ',left,' need ',top0
         call lnkerr('m511: gofish')
      endif
      ioff=0
      call vylm(vlmax,vnomega,ioff,itch(angpts),test,vnomega,nlm,
     $          itch(ylm),aitch(ptlm),itch(ctheta),itch(phi),
     $          itch(cphi),itch(sphi),itch(plm),itch(scr))
      if(debug) then
         do 113 i=1,nlm
            write(iout,*) 'ylm;l,m',i,
     $         (itch(ylm+(i-1)*vnomega+k-1),k=1,vnomega)
  113    continue
      endif
c
c     --- test orthonormality
c         ylmerr is the largest deviation from orthonormality
c         and will be used as a measure of the "machine zero".
      call tstylm(vlmax,vnomega,nlm,itch(ylm),itch(angwt),
     $            aitch(ptlm),itch(scr),ylmerr)
c      write(iout,*) 'maximum ylm error',ylmerr
      if(timeit) then
         call timing(dum4,dum5,dum6)
         timpre=timpre+dum4-dum1
      endif
c
c
c
c     --- generate the potential originating from this charge 
c         distribution.
      if(dogrdio) then
         if (.not.called) then
            call iosys('create real "poisson voronoi weights" on '
     $           //grdfil(1:3),pmxgrd*nat,0,0,' ')
         endif
c         call iosys('rewind "poisson voronoi weights" on '//grdfil(1:3),
c     $              0,0,0,' ')
      endif
c
c     --- prepare density matrix.
      d2=top0
      top0=d2+nnp
      if (top0 .gt. left) then
         write(iout,*)'not enough core for d2!'
         write(iout,*)'have ',left,' need ',top0
         call lnkerr('m511: gofish')
      endif
      call vmove(itch(d2),d,nnp)
      if(zdiag) then
c        --- prepare a density matrix in which the atomic "diagonal" blocks
c            are zeroed.
         do 29 iatom=1,nat
            jatom=iatom
            do 27 itype=1,nbtype
               if (noprim(iatom,itype).le.0) go to 27
               if (iatom.ne.jatom) then
                  jtypmx=nbtype
               else
                  jtypmx=itype
               end if
               do 26 jtype=1,jtypmx
                  if (noprim(jatom,jtype).le.0) go to 26
c
                  nconti=nocont(iatom,itype)
                  ncontj=nocont(jatom,jtype)
                  lenblk=nocart(itype)*nocart(jtype)
                  call rzero(itch(top0),nconti*ncontj*lenblk)
                  call put1el(itch(d2),itch(top0),start,iatom,jatom,
     $                        itype,jtype,nconti,ncontj,nnp,lenblk,nat,
     $                        nbtype,nobf)
   26          continue
   27       continue
   29    continue
   30    continue
         if(debug) then
            write(iout,*) 'nnp',nnp
            write(iout,*) 'd2',(itch(d2+i-1),i=1,nnp)
         endif
      else if(zoffdiag) then
c        --- work only with the diagonal density.
         do 39 iatom=1,nat
            do 38 itype=1,nbtype
               if (noprim(iatom,itype).le.0) go to 38
               do 37 jatom=1,iatom
                  if (iatom.ne.jatom) then
                     jtypmx=nbtype
                  else
c                    note we skip out here.
                     goto 37
                     jtypmx=itype
                  end if
                  do 36 jtype=1,jtypmx
                    if (noprim(jatom,jtype).le.0) go to 36
c
                     nconti=nocont(iatom,itype)
                     ncontj=nocont(jatom,jtype)
                     lenblk=nocart(itype)*nocart(jtype)
                     call rzero(itch(top0),nconti*ncontj*lenblk)
                     call put1el(itch(d2),itch(top0),start,iatom,jatom,
     $                        itype,jtype,nconti,ncontj,nnp,lenblk,nat,
     $                        nbtype,nobf)
   36             continue
   37          continue
   38       continue
   39    continue
   40    continue
         if(debug) then
            write(iout,*) 'nnp',nnp
            write(iout,*) 'd2',(itch(d2+i-1),i=1,nnp)
         endif
      endif
c zero these arrays, because we might be parallel and then have to do a 
c global sum
      call rzero(rpts,vradial*nat)
      call rzero(vlm,vradial*nlm*nat)
      call rzero(y2,vradial*nlm*nat)
c
c     --- solve the poisson equation.
      next=nxtval(nproc)+1
      do 80 iatom=1,nat
         if (iatom.ne.next) goto 80
         if (grdtyp(iatom).eq.'nogrid') goto 80
c        --- put down radial quadrature
         if(timeit) then
            call timing(dum1,dum2,dum3)
         endif
         if(rtox) then
            rwt=top0
            nrings=(vradial/(vncrule-1))+1
            jacob=rwt+nrings*(vncrule*(vncrule-1))
            npring=wpadti(jacob+vradial)
            top=iadtwp(npring+vradial)
            if (top .gt. left) then
               write(iout,*)'not enough core for cotes !'
               write(iout,*)'have ',left,' need ',top
               call lnkerr('m511: gofish')
            endif
            rhomx=rhomax('bragg',ian(iatom))
            call cotes(rhomx,vradial,vncrule,nrings,
     $           rpts(1,iatom),itch(rwt),itch(jacob),
     $           aitch(npring),nrings,nr,nwtot)
c           --- temporary fix, add a large number at end to signal infinity.
c            rpts(vradial+1,iatom)=1.0d20
         else
c           this option is probably broken now that rpts is passed.
c            rwt=top0
c            npring=wpadti(rwt+mxrings*mxcotes*(mxrings*mxcotes-1))
c            scr=iadtwp(npring+mxrings)
c            top=scr+mxrings+1
c            if (top .gt. left) then
c               write(iout,*)'not enough core for newton !'
c               write(iout,*)'have ',left,' need ',top
c               call lnkerr('m511: gofish')
c            endif
c            call newton(ian(iatom),nrings,rpts(1,iatom),itch(rwt),
c     $                  itch(scr),aitch(npring),mxrings,nr,nwtot)
c            write(iout,*)'    newton'
         endif
         if(debug) then
            write(iout,*) 'nrings,nr',nrings,nr
            write(iout,*) 'n per ring',(aitch(npring+k),k=0,nrings-1)
            write(iout,*) 'rpts',(rpts(k,iatom),k=1,nr)
            k=0
            do 114 i=1,nrings
               write(iout,*) 'rwts',(itch(rwt+k),k=0,
     $              aitch(npring+i-1)*(aitch(npring+i-1)-1)-1)
               k=k+aitch(npring+i-1)*(aitch(npring+i-1)-1)  
 114        continue
         endif
c
c        --- combine the radial and angular points into grid
         grid=top
         ngrid=nr*vnomega
         gridwts=grid+3*ngrid
         vptrad=wpadti(gridwts+ngrid)
         top=iadtwp(vptrad+nr+1)
         if (top .gt. left) then
            write(iout,*)'not enough core for catenat !'
            write(iout,*)'have ',left,' need ',top
            call lnkerr('m511: gofish')
         endif
         call catenat(nr,ngrid,vnomega,vnomega,rpts(1,iatom),itch(rwt),
     $        itch(angpts),itch(angwt),itch(grid),
     $        itch(gridwts),aitch(vptrad),c(1,iatom))
c     
c        --- find voronoi wts for this atom
         tvwts=top
c        temporary sleaze
         twts=tvwts+max(pmxgrd,ngrid)
         top=twts+ngrid
         if (top .gt. left) then
            write(iout,*)'not enough core for voronoi !'
            write(iout,*)'have ',left,' need ',top
            call lnkerr('m511: gofish')
         endif
         if(dogrdio.and.didvor(iatom)) then
            call iosys('read real "poisson voronoi weights" from '
     $           //grdfil(1:3),pmxgrd,itch(tvwts),pmxgrd*(iatom-1),' ')
         else
            call voronoi(nat,c,itch(grid),itch(twts),ngrid,ngrid,
     $                   iatom,itch(tvwts),rnuc,amu,
     $                   pwtx,rr,adjust,radii)
            if (dogrdio) then
               call iosys('write real "poisson voronoi weights" to '
     $              //grdfil(1:3),pmxgrd,itch(tvwts),
     $              pmxgrd*(iatom-1),' ')
            endif
            didvor(iatom)=.true.
         endif
         if(debug) then
            write(iout,*) 'vwts',(itch(tvwts+k),k=0,ngrid-1)
         endif
         if(timeit) then
            call timing(dum4,dum5,dum6)
            timagrd=timagrd+dum4-dum1
         endif
c        --- put down the basis functions on the grid.
         if(timeit) then
            call timing(dum1,dum2,dum3)
         endif
         ioff=0
         phibf=top
         grad=phibf
         hess=grad
         dengrad=hess
         scrgrd=phibf+3*ngrid*nbf
         rsq=max(scrgrd+ngrid,scrgrd+minesz*nbf)
         s=rsq+ngrid
         r=s+ngrid*mxcont
         t=r
         xyzpow=r+ngrid*mxcont
         top=xyzpow+ngrid*3*max(bigl,2)
         if (top .gt. left) then
            write(ierr,*)'Node ',mynodeid,
     $           ' Not enough core for bfgrd (need ',top,
     $           ' have ',left,')'
            call plnkerr('oops in gofish',9900)
         endif
         call bfgrd(c,ex,cont,ptprim,noprim,nocont,ptcont,
     $              nat,nprim,ntypes,nbtype,nnp,ncont,
     $              start,nbf,nocart,nobf,maxmom,minmom,mintyp,
     $              nx,ny,nz,itch(grid),ngrid,
     $              itch(xyzpow),itch(scrgrd),itch(rsq),itch(s),
     $              itch(r),itch(t),itch(phibf),itch(grad),itch(hess),
     $              .false.,.false.,mxcont,
     $              ngrid,ngrid,ioff,maxl,bfcut)
         if(timeit) then
            call timing(dum4,dum5,dum6)
            timbf=timbf+dum4-dum1
         endif
c
c        --- form density on the grid for this atom,
c            must worry about ndmat here someday
         call timing(dum1,dum2,dum3)
         dengrid=rsq
         dengrad=dengrid
         phibar=dengrid+ngrid
         gradbar=phibar
         top=gradbar+nbf
         if (top .gt. left) then
            write(iout,*)'not enough core for gridden !'
            write(iout,*)'have ',left,' need ',top
            call lnkerr('m511: gofish')
         endif
         call rzero(itch(dengrid),ngrid)
c        must modify the one below to do d(1,dmat)
         call gridden(nbf,nnp,ngrid,ngrid,itch(d2),itch(phibf),
     $                itch(grad),itch(dengrid),itch(dengrad),
     $                minesz,itch(scrgrd),itch(phibar),itch(gradbar),
     $                dmcut,.false.)
         if(debug) then
            write(iout,*) 'dengrid',(itch(dengrid+k),k=0,ngrid-1)
         endif
         call timing(dum4,dum5,dum6)
         timden=timden+dum4-dum1
         if(timeit) then
            call timing(dum1,dum2,dum3)
         endif
c
c        --- multiply by the voronoi weights to get the atomic
c            component of the density.
         call vmul(itch(dengrid),itch(dengrid),itch(tvwts),ngrid)
         if(debug) then
            write(iout,*) 'weighted',(itch(dengrid+k),k=0,ngrid-1)
         endif
c
c        --- decompose density into ylm components
         rholm=phibar
         scr=rholm+nr*nlm
         top=scr+nr*nlm
         if (top .gt. left) then
            write(iout,*)'not enough core for ftoylm (2) !'
            write(iout,*)'have ',left,' need ',top
            call lnkerr('m511: gofish')
         endif
         call ftoylm(itch(dengrid),itch(rholm),itch(ylm),
     $               nr,vnomega,nlm,itch(angwt),itch(scr),ylmerr)
c         if(debug) then
            write(iout,*) 'rholm'
            do 111 i=1,nlm
               write(iout,*) 'lm',(itch(rholm+(i-1)*nr+k),k=0,nr-1)
  111       continue
c         endif
         if(timeit) then
            call timing(dum4,dum5,dum6)
            tim2ylm=tim2ylm+dum4-dum1
         endif
c
c        --- and compute total charge to see how well the grid does.
c         do 112 i=1,nlm
c            call vmul(scr,itch(rholm+(i-1)*nr),rpts(1,iatom),nr)
c            call vmul(scr,scr,rpts(1,iatom),nr)
c            charge(iatom,1)=charge(iatom,1)
c     $                     +sdot(nr,scr,1,itch(jacob),1)
c  112    continue
c
c        --- solve the radial equations
         if(timeit) then
            call timing(dum1,dum2,dum3)
         endif
         int0f=scr
         int0b=int0f+nlm
         scr=int0b+nlm
         tmp=scr+nr
         j=tmp+nr
         y=j+nr*(vlmax+1)
         top=y+nr*(vlmax+1)
         if (top .gt. left) then
            write(iout,*)'not enough core for rsolver !'
            write(iout,*)'have ',left,' need ',top
            call lnkerr('m511: gofish')
         endif
c
c        --- spline fit the radial solutions
c$$$         ioff=0
c$$$         do 60 lm=1,nlm
c$$$            if(debug) then
c$$$               write(iout,*) 'vlm',(vlm(i,lm,iatom),i=1,nr)
c$$$            endif
c$$$            if(debug) then
c$$$              write(iout,*) 'sum to fit',lm
c$$$              do 61 i=1,nr
c$$$                 write(iout,*) rpts(i,iatom),vlm(i,lm,iatom)
c$$$   61         continue
c$$$            endif
c$$$            call spline3(rpts(1,iatom),vlm(1,lm,iatom),nr,1.0d+31,
c$$$     $                   1.0d+31,y2(1,lm,iatom),itch(scr))
c$$$            ioff=ioff+nr
c$$$   60    continue
c$$$         write(ierr,*)'Here we are, boss node=',mynodeid
c$$$         write(4,*)'Here we are, boss node=',mynodeid
c$$$         if (dumpfit) then
c$$$            do 6060 l=0,vlmax
c$$$               do 6061 m=-l,l
c$$$                  write(4,*)'atom, l,m',iatom,l,m
c$$$                  write(4,*)'sum to fit'
c$$$                  do 6062 i=1,nr
c$$$                     write(4,*) rpts(i,iatom),vlm(i,
c$$$     $                    aitch(ptlm+l+(m+vlmax)*(vlmax+1)),iatom)
c$$$ 6062             continue 
c$$$ 6061          continue 
c$$$ 6060       continue 
c$$$         endif
         call rsolver(nrings,nr,aitch(npring),vlmax,nlm,
     $                aitch(ptlm),rpts(1,iatom),itch(rwt),nwtot,
     $                itch(rholm),itch(j),itch(y),itch(scr),
     $                itch(tmp),vlm(1,1,iatom),itch(int0f),itch(int0b))
         write(iout,*)' vlms'
         do 666 i=1,nlm
            write(iout,*)' channel ',i
            write(iout,*)(vlm(ig,i,iatom),ig=1,nr)
 666     continue 
         if(timeit) then
            call timing(dum4,dum5,dum6)
            timrsolv=timrsolv+dum4-dum1
         endif
         call timing(dum1,dum2,dum3)
c
c new idea:
c  Instead of doing all atoms, spline fitting, mapping onto every grid, 
c and then repartitioning into atomic contributions a la Becke, why not
c do a straight cotes/lebedev quadrature on the V(n) that we already have?
c
c  We already have the basis on this atom's grid.  We have this atom's radial
c  quadrature points, and the ulm's evaluated exactly at them.  We have
c  the ylm's on the angular grid.  Just do a straight sum over grid points
c  with the weights, BUT LEAVE OUT THE VORONOI WEIGHTS.  The J matrix will
c  then be a simple sum over all of these contributions, which global sum is 
c  done in the calling routine.
c
c        ---form total v(n) on this atom's grid
         v=int0f
         top=v+vnomega*vradial
         call ftogrid0(itch(ylm),aitch(ptlm),vnomega,vlmax,
     $        vradial,rpts(1,iatom),vlm(1,1,iatom),nlm,itch(v))
         call timing(dum4,dum5,dum6)
         tim2grd=tim2grd+dum4-dum1
c$$$         write(iout,*)' v is:'
c$$$         do 1010 ig=1,ngrid
c$$$            write(iout,*) 'v(',ig,')=',itch(ig+v-1)
c$$$ 1010    continue 
c        ---form jmatrix integrals using basis on this grid and add to
c           current j matrix
         call timing(dum1,dum2,dum3)
c
c the cotes rwts are *NOT* integration weights.  It turns out that if
c we take the cotes *jacob* array, divide it by nr, and multiply it by
c pts**2 we'll get back the euler-maclaurin weights.  Do that here
c
         do 1111 ig=0,nr-1
            itch(rwt+ig)=itch(jacob+ig)*rpts(ig+1,iatom)**2/float(nr)
c            if (mynodeid.eq.0)write(iout,*)'rwt(',ig+1,')=',
c     $           itch(rwt+ig)
 1111     continue 
         call catenat(nr,ngrid,vnomega,vnomega,rpts(1,iatom),itch(rwt),
     $        itch(angpts),itch(angwt),itch(grid),
     $        itch(gridwts),aitch(vptrad),c(1,iatom))
c$$$         do 1112 ig=1,ngrid
c$$$            write(iout,*)'wt(',ig,')=',itch(gridwts+ig-1)
c$$$ 1112    continue 
c
c sleaze because fmkclos actually integrates fout(.,2) not fout(.,1)
c
         fout=v-ngrid
         tea=top
         scr=tea+minesz*nbf
         queue=scr+minesz
         phibar=queue+3*minesz
         gradbar=phibar
         top=gradbar+nbf
         call fmkclos(ngrid,ngrid,minesz,nbf,nnp,.false.,.false.,
     $        itch(gridwts),itch(scr),itch(phibar),itch(gradbar),
     $        itch(queue),itch(tea),itch(fout),itch(phibf),
     $        itch(grad),itch(dengrad),jmat,kmcut)
         write(iout,*)' after atom ',iatom,' integration, Jmat is'
         call print(jmat,nnp,nbf,iout)
         call timing(dum4,dum5,dum6)
         timfmk=timfmk+dum4-dum1
         next=nxtval(nproc)+1
   80 continue
      next=nxtval(-nproc)
      if(timeit) then
         write(iout,*) 'timpre,timagrd,timbf,timden',
     $                  timpre,timagrd,timbf,timden
         write(iout,*) 'tim2ylm,timrsolv,timatg,tim2grd,timfmk',
     $                  tim2ylm,timrsolv,timatg,tim2grd,timfmk
      endif
      write(4,*) '      node ',mynodeid,' gofish gridden', timden
      write(4,*) '      node ',mynodeid,' gofish v to grid', tim2grd
      write(4,*) '      node ',mynodeid,' gofish form J', timfmk

c
c     --- make sure we know we've been here.
      called=.true.
c
c
      return
      end
