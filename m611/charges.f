*deck @(#)charges.f	5.4 11/28/95
      subroutine charges(itch,d,nbf,nnp,ndmat,nat,nradial,lmax,
     $                   order,dmcut,dencut,adjust,aitch,
     $                   ian,c,ex,cont,ptprim,noprim,nocont,ptcont,
     $                   mxcont,nprim,ntypes,nbtype,ncont,
     $                   start,nocart,nobf,maxmom,minmom,mintyp,
     $                   nx,ny,nz,lenxyz,charge,size,maxl,bigl,minesz,
     $                   ops)
c***begin prologue     charges.f
c***date written       940304      (yymmdd)  
c***revision date      11/28/95      
c   july 1, 1994
c      modifying gofish.f to do atomic charges
c
c***keywords           charges, ylm, density, atomic size
c***author             martin, richard(lanl) 
c***source             @(#)charges.f	5.4   11/28/95
c***purpose            decomposes the density into atomic ylm components
c***description        
c
c***references
c
c***routines called
c
c***end prologue       charges.f
      implicit none
c     --- input variables ---
      integer nbf,nnp,ndmat,nat
      integer nprim,ntypes,nbtype,ncont,mxcont
      integer bigl,lenxyz,minesz
      integer nradial,lmax,order
      logical adjust
      character*(*) ops
      real*8 dencut,dmcut
c     --- input arrays (unmodified) ---
      integer ian(nat)
      integer ptprim(nat,ntypes),noprim(nat,ntypes),nocont(nat,ntypes)
      integer ptcont(nat,ntypes),start(nat,ntypes),nocart(0:*)
      integer nobf(ntypes),maxmom(ntypes),minmom(ntypes),mintyp(ntypes)
      integer nx(lenxyz),ny(lenxyz),nz(lenxyz)
      real*8 d(nnp,ndmat)
      real*8 c(3,nat),ex(nprim),cont(ncont)
c     --- input arrays (scratch) ---
      real*8 itch(*)
      integer aitch(*)
      integer maxl(nat)
c     --- output arrays ---
      real*8 charge(nat,ndmat),size(nat)
c     --- output variables ---
c     --- scratch arrays ---
c     --- local variables ---
      integer iadtwp,wpadti
      integer nomega,ngrid,lebord
      integer rmax,mxgrd,gradwt,radshls
      integer radii,top,top0,grid,ptrad,akl
      integer phibf,scrgrd,rsq,r,s,t,grad,hess,xyzpow
      integer dengrid,dengrad
      integer vwts,wts,rnuc,amu,pwtx,rr,pr
      integer dofr,y2,u,scr
      integer iatom,ioff
      integer k,dmat
      integer inp,iout
      integer phibar,gradbar
      logical debug
      real*8 zero
      real*8 chg
      real*8 sdot,atan
      real*8 one,two,four,pi
c
      parameter (debug=.false.)
      parameter (zero=0.0d+00,two=2.0d+00,four=4.0d+00)
      data one/1.0d0/
c
      common/io/inp,iout
c      
      integer angsiz
c
c     --- initialize the output array
      call rzero(charge,nat*ndmat)
      pi=four*atan(one)
      top0=1
c
c     --- generate the atomic component corresponding to this charge 
c         distribution.
c         the standard grids take at most 60 radial points.
      rmax=nradial
      lebord=2*lmax+1
      nomega=angsiz(lebord)
      mxgrd=max(nradial*nomega,4000)
      do 100 iatom=1,nat
c
c        --- generate the grid and weights for this atom
         grid=top0
         wts=grid+mxgrd*3
         gradwt=wts
         ptrad=wpadti(wts+mxgrd)
         rnuc=iadtwp(ptrad+rmax)
         amu=rnuc+nat*nat
         radii=amu+nat*nat
         akl=radii+nat
         pwtx=akl+nat*nat
         rr=pwtx+nat
         vwts=rr+nat
         top=vwts+mxgrd
         call mkatmg(c,ian,itch(grid),itch(wts),rmax,lebord,nomega,
     $               nradial,nat,ngrid,mxgrd,itch(vwts),itch(rnuc),
     $               itch(amu),itch(pwtx),itch(rr),adjust,itch(radii),
     $               itch(akl),'general',radshls,aitch(ptrad),
     $               .false.,itch(gradwt),iatom)
         if(debug) then
            write(iout,*) 'wts',(itch(wts+k),k=0,ngrid-1)
         endif
c
c        --- put down the basis functions on the grid.
         phibf=rnuc
         grad=phibf
         hess=grad
         scrgrd=phibf+3*ngrid*nbf
c        use scrgrd for scratch in gridden too.
         rsq=max(scrgrd+ngrid,scrgrd+minesz*nbf)
         s=rsq+ngrid
         r=s+ngrid*mxcont
         t=r
         xyzpow=r+ngrid*mxcont
         top=xyzpow+ngrid*3*max(bigl,2)
         ioff=0
c        assume there is enough left
         call bfgrd(c,ex,cont,ptprim,noprim,nocont,ptcont,
     $              nat,nprim,ntypes,nbtype,nnp,ncont,
     $              start,nbf,nocart,nobf,maxmom,minmom,mintyp,
     $              nx,ny,nz,itch(grid),mxgrd,
     $              itch(xyzpow),itch(scrgrd),itch(rsq),itch(s),
     $              itch(r),itch(t),itch(phibf),itch(grad),itch(hess),
     $              .false.,.false.,mxcont,ngrid,ngrid,ioff,maxl)
c
c        --- form density on the grid for this atom,
         dengrid=rsq
         dengrad=rsq
         phibar=dengrid+ndmat*ngrid
         gradbar=phibar
         top=gradbar+minesz*nbf
         do 10 dmat=1,ndmat
            call rzero(itch(dengrid),ngrid)
            call gridden(nbf,nnp,ngrid,ngrid,d(1,dmat),itch(phibf),
     $                   itch(grad),itch(dengrid),itch(dengrad),
     $                   minesz,itch(scrgrd),itch(phibar),
     $                   itch(gradbar),dmcut,.false.)
            if(debug) then
               write(iout,*) 'dengrid',(itch(dengrid+k),k=0,ngrid-1)
            endif
            chg=sdot(ngrid,itch(dengrid),1,itch(wts),1)
            if(dmat.eq.1) then
               charge(iatom,dmat)=two*chg
            else
               charge(iatom,dmat)=chg
            endif
c           increment dengrid for open shell density
            dengrid=dengrid+ngrid
   10    continue
c 
c        --- determine an estimate of the 'atomic' size.
c            this is done by integrating the density over the unit sphere
c            for each radial shell, then spline fitting the result.
c            reset dengrid to closed shell starting address
         dengrid=rsq
         dofr=dengrid+ndmat*ngrid
         pr=dofr+radshls+1
         y2=pr+radshls+1
         u=y2+radshls+1
         scr=y2+radshls+1
         top=scr+ngrid 
         call sizes(ngrid,radshls,itch(dengrid),itch(wts),aitch(ptrad),
     $              ndmat,itch(pr),itch(dofr),itch(y2),itch(u),
     $              itch(scr),itch(grid),mxgrd,
     $              charge(iatom,1)+charge(iatom,2),size(iatom))
  100 continue
c
c
      return
      end
