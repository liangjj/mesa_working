*deck @(#)gofish2.f	5.1  11/28/95
      subroutine gofish2(itch,d,nbf,nnp,jmat,ncoul,
     $     ndmat,nat,mxgrd,
     $     dmcut,dencut,
     $     mxgbsiz,aitch,ian,
     $     c,ex,cont,ptprim,noprim,nocont,ptcont,
     $     mxcont,nprim,ntypes,nbtype,ncont,
     $     start,nocart,nobf,maxmom,minmom,mintyp,
     $     nx,ny,nz,xyzgrid,grdwts,charge,maxl,bigl,ops,
     $     left,vwts,rnuc,amu,pwtx,rr,radii,akl,ptrad,rmax,lmax,
     $     nomega,nradial,grdtyp,adjust,minesz,vlmax,vradial,
     $     vncrule,nlm,rpts,vlm,y2)
c***begin prologue     gofish2.f
c***date written       940304      (yymmdd)  
c***revision date      4/25/95      
c
c***keywords           poisson, coulomb, potential, density
c***author             martin, richard(lanl) 
c***source             @(#)gofish2.f	5.1   11/28/95
c***purpose            generates the coulomb potential from the density
c***description        
c                      solves poisson equation for v, given rho
c                         (del**2) v = rho
c     
c                      this is accomplished by projecting the total
c                      density into single-center pieces and a remainder.
c                      the remainder is projected onto the atomic grids of 
c                      Becke, thereby reducing the problem to a series of 
c                      atomic poisson problems.  
c
c                      for each atom, the density is decomposed into
c                      spherical harmonic components and a radial equation
c                      is solved for the potential originating from that
c                      component. the radial equation is converted into an
c                      integral equation using the appropriate green's 
c                      function, and solved via newton-core quadrature.
c                      
c                      it should be noted that the grid used to solve the
c                      atomic problem may be  different from that used 
c                      to represent the resulting potential.
c
c***references
c
c***routines called
c
c***end prologue       gofish2.f
      implicit none
c     --- input variables ---
      integer nbf,nnp,ncoul,ndmat,nat,mxgrd,mxgbsiz
      integer nprim,ntypes,nbtype,ncont,mxcont
      integer rmax,lmax,nradial,bigl
      integer left,minesz,nlm
      integer vlmax,vradial,vncrule
      real*8 dencut,dmcut
c     --- input arrays (unmodified) ---
      character*(*) ops
      character*(*) grdtyp(nat)
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
c     --- local variables ---
      logical adjust,rtox,dogrdio
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
      integer i,ioff,joff,iatom,jatom,ngj
      integer lm
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
      logical debug,oddblock
      logical timeit
      real*8 ylmerr,zero,one
      real*8 test(3)
      real*8 dum1,dum2,dum3,dum4,dum5,dum6
      real*8 timpre,timagrd,timbf,timden,tim2ylm,timrsolv
      real*8 timatg,tim2grd,timfmk
      real*8 sdot
c
      parameter (debug=.false.,timeit=.true.)
      parameter (rtox=.true.,dogrdio=.false.)
      parameter (mxrings=100,mxcotes=9)
      parameter (zero=0.0d+00,one=1.0d+00)
      data timpre,timagrd,timbf,timden,tim2ylm,timrsolv/6*0.0d0/
      data timatg,tim2grd,timfmk/3*0.0d0/
c
      common/io/inp,iout
c      
      integer angsiz
      real*8 rhomax,rhomx 
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
         call lnkerr('m515: gofish')
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
         call lnkerr('m515: gofish')
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
      do 80 iatom=1,nat
c
c        --- put down radial quadrature
            if(timeit) then
               call timing(dum1,dum2,dum3)
            endif
         if(rtox) then
            rwt=ctheta
            nrings=(vradial/(vncrule-1))+1
            jacob=rwt+nrings*(vncrule*(vncrule-1))
            npring=wpadti(jacob+vradial)
            top=iadtwp(npring+vradial)
            if (top .gt. left) then
               write(iout,*)'not enough core for cotes !'
               write(iout,*)'have ',left,' need ',top
               call lnkerr('m515: gofish')
            endif
            rhomx=rhomax('bragg',ian(iatom))
            call cotes(rhomx,vradial,vncrule,nrings,
     $                 rpts(1,iatom),itch(rwt),itch(jacob),
     $                 aitch(npring),nrings,nr,nwtot)
         else
c           this option is probably broken now that rpts is passed.
c            rwt=ctheta
c            npring=wpadti(rwt+mxrings*mxcotes*(mxrings*mxcotes-1))
c            scr=iadtwp(npring+mxrings)
c            top=scr+mxrings+1
c            if (top .gt. left) then
c               write(iout,*)'not enough core for newton !'
c               write(iout,*)'have ',left,' need ',top
c               call lnkerr('m515: gofish')
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
     $             aitch(npring+i-1)*(aitch(npring+i-1)-1)-1)
               k=k+aitch(npring+i-1)*(aitch(npring+i-1)-1)  
  114       continue
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
            call lnkerr('m515: gofish')
         endif
         call catenat(nr,ngrid,vnomega,vnomega,rpts(1,iatom),itch(rwt),
     $               itch(angpts),itch(angwt),itch(grid),
     $               itch(gridwts),aitch(vptrad),c(1,iatom))
c
c        --- compute the atomic block analytically.
c            need an atomic spline representation of l=0.
c       STILL A PROBLEM
c
c        --- find voronoi wts for this atom
         tvwts=top
         twts=tvwts+ngrid
         top=twts+ngrid
         if (top .gt. left) then
            write(iout,*)'not enough core for voronoi !'
            write(iout,*)'have ',left,' need ',top
            call lnkerr('m515: gofish')
         endif
         call voronoi(nat,c,itch(grid),itch(twts),ngrid,ngrid,
     $                iatom,itch(tvwts),rnuc,amu,
     $                pwtx,rr,adjust,radii)
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
            write(iout,*)' Not enough core for bfgrd (need ',top,
     $           ' have ',left,')'
            call lnkerr('oops in gofish')
         endif
         call bfgrd(c,ex,cont,ptprim,noprim,nocont,ptcont,
     $              nat,nprim,ntypes,nbtype,nnp,ncont,
     $              start,nbf,nocart,nobf,maxmom,minmom,mintyp,
     $              nx,ny,nz,itch(grid),ngrid,
     $              itch(xyzpow),itch(scrgrd),itch(rsq),itch(s),
     $              itch(r),itch(t),itch(phibf),itch(grad),itch(hess),
     $              .false.,.false.,mxcont,
     $              ngrid,ngrid,ioff,maxl)
            if(timeit) then
               call timing(dum4,dum5,dum6)
               timbf=timbf+dum4-dum1
            endif
c
c        --- form density on the grid for this atom,
c            must worry about ndmat here someday
            if(timeit) then
               call timing(dum1,dum2,dum3)
            endif
         dengrid=rsq
         dengrad=dengrid
         phibar=dengrid+ngrid
         gradbar=phibar
         d2=gradbar+nbf
         top=d2+nnp
         if (top .gt. left) then
            write(iout,*)'not enough core for gridden !'
            write(iout,*)'have ',left,' need ',top
            call lnkerr('m515: gofish')
         endif
         call rzero(itch(dengrid),ngrid)
         call vmove(itch(d2),d(1,1),nnp)
c        must modify the one below to do d(1,dmat)
c        modify the density matrix to omit the diagonal blocks.
         call gridden(nbf,nnp,ngrid,ngrid,itch(d2),itch(phibf),
     $                itch(grad),itch(dengrid),itch(dengrad),
     $                minesz,itch(scrgrd),itch(phibar),itch(gradbar),
     $                dmcut,.false.)
         if(debug) then
            write(iout,*) 'dengrid',(itch(dengrid+k),k=0,ngrid-1)
         endif
            if(timeit) then
               call timing(dum4,dum5,dum6)
               timden=timden+dum4-dum1
            endif
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
            call lnkerr('m515: gofish')
         endif
         call ftoylm(itch(dengrid),itch(rholm),itch(ylm),
     $               nr,vnomega,nlm,itch(angwt),itch(scr),ylmerr)
         if(debug) then
            write(iout,*) 'flm'
            do 111 i=1,nlm
               write(iout,*) 'lm',(itch(rholm+(i-1)*nr+k),k=0,nr-1)
  111       continue
         endif
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
            call lnkerr('m515: gofish')
         endif
         call rsolver(nrings,nr,aitch(npring),vlmax,nlm,
     $                aitch(ptlm),rpts(1,iatom),itch(rwt),nwtot,
     $                itch(rholm),itch(j),itch(y),itch(scr),
     $                itch(tmp),vlm(1,1,iatom),itch(int0f),itch(int0b))
         if(debug) then
            write(iout,*) 'iatom,int0f',iatom,(itch(int0f+i-1),i=1,nlm)
            write(iout,*) 'iatom,int0b',iatom,(itch(int0b+i-1),i=1,nlm)
            ioff=0
c            do 71 lm=1,nlm
c               call vadd(itch(vlm+ioff),
c     $                   itch(vlm+ioff),itch(falm+ioff),nr)
c               write(iout,*) 'ulm diffs',lm
c                  do 72 i=1,nr
c                     diff=abs(itch(falm+ioff+i-1)-itch(vlm+ioff+i-1))
c                     if(diff.gt.1.0d-08) then
c                       write(iout,*) 'pt,diff',i,diff
c                     endif
c   72             continue
c               ioff=ioff+nr
c   71       continue
         endif
c 
c        --- spline fit the radial solutions
c        --- don't forget to combine the analytical and numerical radial
c            contributions
c        will want to modify this to include atomic density.
         ioff=0
         do 60 lm=1,nlm
            if(debug) then
               write(iout,*) 'vlm',(vlm(i,lm,iatom),i=1,nr)
            endif
c            call vadd(itch(vlm+ioff),itch(vlm+ioff),itch(falm+ioff),nr)
ctmp            call vmove(itch(vlm+ioff),itch(falm+ioff),nr)
             if(debug) then
               write(iout,*) 'sum to fit',lm
               do 61 i=1,nr
                  write(iout,*) rpts(i,iatom),vlm(i,lm,iatom)
   61          continue
             endif
             call spline3(rpts(1,iatom),vlm(1,lm,iatom),nr,1.0d+31,
     $                    1.0d+31,y2(1,lm,iatom),itch(scr))
            ioff=ioff+nr
   60    continue
            if(timeit) then
               call timing(dum4,dum5,dum6)
               timrsolv=timrsolv+dum4-dum1
            endif
   80 continue
c
c     --- we now have the potential decomposed into atomic contributions
c         and fit to a cubic spline: rpts,vlm,y2.
c
c        --- put the solution back onto the full grid
c            note that here we use the original grid.
      top0=1
      if(dogrdio) then
         grdsiz=1
         top0=iadtwp(grdsiz+nat)
         call iosys('read "external atomic grid size" from rwf',
     $              nat,aitch(grdsiz),0,' ') 
         call iosys('rewind "external grid" on rwf',0,0,0,' ')
         call iosys('rewind "external grid weights" on rwf',
     $               0,0,0,' ')
      endif
      do 100 jatom=1,nat
c
c        --- generate grid, calc number of grid blocks
         if(timeit) then
            call timing(dum1,dum2,dum3)
         endif
         if(dogrdio) then
            jngrid=aitch(grdsiz+jatom-1)
            call iosys('read real "external grid" from rwf'
     $                 //' without rewinding',mxgrd*3,xyzgrid,0,' ')
            call iosys('read real "external grid weights" from rwf'
     $                 //' without rewinding',mxgrd,grdwts,0,' ')
         else
            call mkatmg(c,ian,xyzgrid,grdwts,rmax,lmax,nomega,
     $                  nradial,
     $                  nat,jngrid,mxgrd,vwts,rnuc,amu,pwtx,rr,adjust,
     $                  radii,akl,grdtyp(jatom),
     $                  radshls,ptrad,.false.,vwts,jatom)
         endif
         if(timeit) then
            call timing(dum4,dum5,dum6)
            timatg=timatg+dum4-dum1
         endif
         ngb=jngrid/mxgbsiz
         oddblock=mod(jngrid,mxgbsiz)
         if (oddblock.ne.0) ngb=ngb+1
         joff=0
         do 90 jgblk=1,ngb
            if (jgblk .ne. ngb .or. oddblock.eq.0 ) then
               ngj=mxgbsiz
            else
               ngj=oddblock
            endif
            v=top0
            fr=v+ngj
            ylm=fr+ngj*nlm
            ptlm=wpadti(ylm+ngj*nlm)
            ctheta=iadtwp(ptlm+(vlmax+1)*(2*vlmax+1))
            phi=ctheta+ngj
            cphi=phi+ngj
            sphi=cphi+ngj
            plm=sphi+ngj
            scr=plm+ngj*(vlmax+1)
            ind=wpadti(scr+ngj)
            top=iadtwp(ind+ngj)
            if (top .gt. left) then
               write(iout,*)'not enough core for ftogrid !'
               write(iout,*)'have ',left,' need ',top
               call lnkerr('m515: gofish')
            endif
c
            call rzero(itch(v),ngj)
            if(timeit) then
               call timing(dum1,dum2,dum3)
            endif
            do 75 iatom=1,nat
               call ftogrid(ngj,vlmax,nlm,mxgrd,joff,
     $                      xyzgrid,c(1,iatom),itch(v),
     $                      itch(fr),itch(ylm),aitch(ptlm),itch(plm),
     $                      itch(ctheta),itch(phi),itch(cphi),
     $                      itch(sphi),itch(scr),nr,rpts(1,iatom),
     $                      vlm(1,1,iatom),y2(1,1,iatom),aitch(ind))
               if(debug) then
                  write(iout,*) 'fitted v'
                  write(iout,*) (itch(v+i-1),i=1,ngj)
               endif
   75       continue
            if(timeit) then
               call timing(dum4,dum5,dum6)
               tim2grd=tim2grd+dum4-dum1
            endif
c
c           --- compute basis functions on this grid block.
            phibf=fr
            grad=phibf
            hess=grad
            dengrad=hess
            scrgrd=phibf+3*mxgbsiz*nbf
            rsq=max(scrgrd+mxgbsiz,scrgrd+minesz*nbf)
            s=rsq+mxgbsiz
            r=s+mxgbsiz*mxcont
            t=r
            xyzpow=r+mxgbsiz*mxcont
            top=xyzpow+mxgbsiz*3*max(bigl,2)
            if(timeit) then
               call timing(dum1,dum2,dum3)
            endif
            call bfgrd(c,ex,cont,ptprim,noprim,nocont,ptcont,
     $                 nat,nprim,ntypes,nbtype,nnp,ncont,
     $                 start,nbf,nocart,nobf,maxmom,minmom,mintyp,
     $                 nx,ny,nz,xyzgrid,mxgrd,
     $                 itch(xyzpow),itch(scrgrd),itch(rsq),itch(s),
     $                 itch(r),itch(t),itch(phibf),itch(grad),
     $                 itch(hess),.false.,.false.,mxcont,
     $                 mxgbsiz,ngj,joff,maxl)
            if(timeit) then
               call timing(dum4,dum5,dum6)
               timbf=timbf+dum4-dum1
            endif
c     NOTE TIS SLEAZE -- CHANGE FOUT IN THE CALL.
            fout=v-mxgbsiz
            tea=scrgrd
            scr=tea+minesz*nbf
            queue=scr+minesz
            phibar=queue+3*minesz
            gradbar=phibar
            top=gradbar+nbf
            if(timeit) then
               call timing(dum1,dum2,dum3)
            endif
            call fmkclos(ngj,mxgbsiz,minesz,nbf,nnp,
     $                   .false.,.false.,
     $                   grdwts,itch(scr),itch(phibar),itch(gradbar),
     $                   itch(queue),itch(tea),itch(fout),
     $                   itch(phibf),itch(grad),itch(dengrad),
     $                   jmat)
            joff=joff+ngj
            if(timeit) then
               call timing(dum4,dum5,dum6)
               timfmk=timfmk+dum4-dum1
            endif
   90    continue
  100 continue
c
c     --- write out the charges to see how good the grid is
c      do 110 iatom=1,nat
c         write(iout,*) 'atom,charge',iatom,charge(iatom,1)
c  110 continue
c     NOTE THAT WE COULD DELAY THE LAST LOOP UNTIL THE K MATRIX.
c
      if(timeit) then
         write(iout,*) 'timpre,timagrd,timbf,timden',
     $                  timpre,timagrd,timbf,timden
         write(iout,*) 'tim2ylm,timrsolv,timatg,tim2grd,timfmk',
     $                  tim2ylm,timrsolv,timatg,tim2grd,timfmk
      endif
c
      return
      end
