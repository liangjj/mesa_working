*deck %W% %G%
      subroutine gofish(itch,d,nbf,nnp,jmat,ncoul,
     $     ndmat,nat,mxgrd,
     $     dmcut,dencut,ngb,
     $     gblksz,mxgblk,mxgbsiz,aitch,ian,
     $     c,ex,cont,ptprim,noprim,nocont,ptcont,
     $     mxcont,nprim,ntypes,nbtype,ncont,
     $     start,nocart,nobf,maxmom,minmom,mintyp,
     $     nx,ny,nz,xyzgrid,charge,maxl,dograd,bigl,ops)
c***begin prologue     %M%
c***date written       940304      (yymmdd)  
c***revision date      %G%      
c
c***keywords           poisson, coulomb, potential, density
c***author             martin, richard(lanl) 
c***source             %W%   %G%
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
c***end prologue       %M%
      implicit none
c     --- input variables ---
      integer nbf,nnp,ncoul,ndmat,nat,mxgrd,mxgbsiz,mxgblk
      integer nprim,ntypes,nbtype,ncont,mxcont
      integer bigl
      real*8 dencut,dmcut
      logical dograd
c     --- input arrays (unmodified) ---
      character*(*) ops
      integer ngb(nat),gblksz(mxgblk,nat)
      integer ian(nat)
      integer ptprim(nat,ntypes),noprim(nat,ntypes),nocont(nat,ntypes)
      integer ptcont(nat,ntypes),start(nat,ntypes),nocart(0:*)
      integer nobf(ntypes),maxmom(ntypes),minmom(ntypes),mintyp(ntypes)
      integer nx(*),ny(*),nz(*)
      real*8 d(nnp,ndmat)
      real*8 c(3,nat),ex(nprim),cont(ncont)
      real*8 xyzgrid(mxgrd,3,nat)
c     --- input arrays (scratch) ---
      real*8 itch(*)
      integer aitch(*)
      integer maxl(nat)
c     --- output arrays ---
      real*8 jmat(nnp,ncoul)
      real*8 charge(nat,ndmat)
c     --- output variables ---
c     --- scratch arrays ---
c     --- local variables ---
      logical adjust,rtox
      integer intkey,ncrule
      integer wpadti,iadtwp
      integer nomega,nr,ngrid,nlm,lmax,lebord
      integer v,radii,angpts,angwt,top,top0,rpts,rwt,grid,gridwts,ptrad
      integer phibf,scrgrd,rsq,r,s,grad,xyzpow
      integer dengrid,dengrad
      integer vwts,wts,rnuc,amu,pwtx,rr
      integer ylm,ptlm,ctheta,phi,cphi,sphi,plm,scr,tmp
      integer rholm
      integer sc,vlm,coef,break,ind,fr
      integer i,ioff,joff,iatom,jatom,ngi,ngj,igblk,jgblk
      integer order,lm,nbreak
      integer coord,pt,ptv
      integer inp,iout
      integer k,jacob
      integer nrings
      integer npring
      integer mxrings,mxcotes
      integer nwtot
      integer j,y
      integer scr1,vnuc,vel,zan,maxcor,nugrid
      integer poten
      integer va,valm,falm,d2
      integer itype,jtype,jtypmx,nconti,ncontj,lenblk
      integer nradial,int0f,int0b,ntest
      integer ptx,pty,ptz
      real*8 r2
      logical debug,compare
      real*8 sdot
      real*8 ylmerr,diff,mxerr,terr,zero
      real*8 test(3)
c
      parameter (debug=.false.,compare=.true.)
      parameter (rtox=.true.)
      parameter (mxrings=100,mxcotes=9)
      parameter (zero=0.0d+00)
c
      common/io/inp,iout
c      
      integer lebpts(29)
      data lebpts/2*0,6,5*0,38,0,50,0,74,0,86,0,110,0,146,0,0,0,194,
     $            6*0/
      real*8 alpha(36)
c                h-he
      data alpha/1.00,0.59,
c                li-ne
     $           3.08,2.06, 1.55,1.23,1.04,0.89,0.77,0.68,
c                na-ar
     $           4.10,3.17, 2.59,2.17,1.89,1.66,1.47,1.33,
c                k,ca
     $           6.27,4.84, 
c                sc-zn
     $           4.59,4.38,4.20,4.01,3.82,3.69,3.53,3.40,3.27,3.16,
c                ga-kr
     $                      2.76,2.44,2.19,1.98,1.81,1.66/
c
 1000 format(5x,'poisson:',
     $      /8x,'maximum l:        ',i3,
     $      /8x,'nradial:          ',i3,
     $      /8x,'newton-cotes rule:',i3,
     $      /8x,'spline-order:     ',i3)
c
c     --- get the parameters to be used for the atomic decompositions
      order=intkey(ops,'poisson=spline-order',3,' ')
      order=order+1
      lmax=intkey(ops,'poisson=lmax',4,' ')
      nradial=intkey(ops,'poisson=nradial',51,' ')
      ncrule=intkey(ops,'poisson=rule',5,' ')
      write(iout,1000) lmax,nradial,ncrule,order-1
c
c     --- make room for the numerical potential on the input grid.
c         allocate core into itch
      zan=1
      v=zan+nat
      top0=v+mxgrd*nat
      call rzero(itch(v),mxgrd*nat)
      call iosys('read real "nuclear charges" from rwf',
     $            -1,itch(zan),0,' ')
c
c     --- if we are going to compute the potential on the grid
c         analytically to compare with the numerical solution,
c         do it now
      if(compare) then
c        calculate the potential analytically
         poten=top0
         top0=poten+mxgrd*nat
         do 11 iatom=1,nat
            ioff=0
            do 10 igblk=1,ngb(iatom)
               ngi=gblksz(igblk,iatom)
               vel=top0
               vnuc=vel+nnp*ngi 
               nugrid=vnuc+nat
               scr=nugrid+3*ngi
               call getscm(0,itch(scr),maxcor,'gofish',0)
c              the routine v0 expects the grid to be (3,ngi)
c              at this point we have grid(mxgrd,3)
c              turn it around.
               k=0
               do 31 pt=1,ngi
                  do 21 coord=1,3
                     itch(nugrid+k)=
     $                  xyzgrid(pt+ioff,coord,iatom)
                     k=k+1
   21             continue
   31          continue
               call v0(c,ex,itch(scr),cont,ptprim,noprim,nocont,ptcont,
     $                 nat,nprim,iadtwp(maxcor),ntypes,nbtype,nnp,ncont,
     $                 start,nbf,itch(zan),nocart,nobf,maxmom,mintyp,
     $                 nx,ny,nz,minmom,' ',itch(vel),itch(vnuc),.false.,
     $                 ' ',ngi,itch(nugrid))
               scr=vnuc
               scr1=scr+nbf*nbf
               joff=0
               do 121 i=1,ngi
                  call trtosq(itch(scr),d,nbf,nnp)
                  call trtosq(itch(scr1),itch(vel+joff),nbf,nnp)
                  itch(poten+ioff+i-1+(iatom-1)*mxgrd)=
     $               sdot(nbf*nbf,itch(scr),1,itch(scr1),1)
                  joff=joff+nnp
  121          continue
               if(debug) then
                  write(iout,*) 'vanal',
     $               (itch(poten+ioff+i+(iatom-1)*mxgrd),i=0,ngi-1)
               endif
               ioff=ioff+ngi
   10       continue
   11    continue
      endif
c
c     --- prepare a density matrix in which the atomic "diagonal" blocks
c         are zeroed.
      d2=top0
      top0=d2+nnp
      call vmove(itch(d2),d,nnp)
c     go to 30
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
               call put1el(itch(d2),itch(top0),start,iatom,jatom,itype,
     $                     jtype,nconti,ncontj,nnp,lenblk,nat,nbtype,
     $                     nobf)
   26       continue
   27    continue
   29 continue
   30 continue
      if(debug) then
         write(iout,*) 'd2',(itch(d2+i-1),i=1,nnp)
      endif
c
c     --- if we are using adjusted cell radii for the voronoi
c         weights on the poisson grid, do it now
      radii=top0
      adjust=.false.
      if(adjust) then
         top0=radii+nat
         do 310 i=1,nat
            itch(radii+i-1)=alpha(ian(i))
  310    continue
      endif
      call rzero(test,3)
c
c     --- generate the potential originating from this charge 
c         distribution.
      do 100 iatom=1,nat
c
c        --- generate lebedev angular points and weights on the 
c            unit sphere.
         lebord=2*lmax+1
         nomega=lebpts(lebord)
         nlm=(lmax+1)*(lmax+1)
c        --- allocate some core
         angpts=top0
         angwt=angpts+3*nomega
         top=angwt+3*nomega
         call lebedev(nomega,itch(angpts),itch(angwt),lebord,nomega)
         if(debug) then
            write(iout,*) 'angtps',(itch(angpts+k),k=0,3*nomega-1)
            write(iout,*) 'angwt',(itch(angwt+k),k=0,nomega-1)
         endif
c
c        --- generate ylm on the lebedev grid
         ylm=top
         ptlm=wpadti(ylm+nomega*nlm)
         ctheta=iadtwp(ptlm+(lmax+1)*(2*lmax+1))
         phi=ctheta+nomega
         cphi=phi+nomega
         sphi=cphi+nomega
         plm=sphi+nomega
         scr=plm+nomega*(lmax+1)
         top=scr+nomega
         ioff=0
         call vylm(lmax,nomega,ioff,itch(angpts),test,nomega,nlm,
     $             itch(ylm),aitch(ptlm),itch(ctheta),itch(phi),
     $             itch(cphi),itch(sphi),itch(plm),itch(scr))
         if(debug) then
            do 113 i=1,nlm
               write(iout,*) 'ylm;l,m',i,
     $            (itch(ylm+(i-1)*nomega+k-1),k=1,nomega)
  113       continue
         endif
c
c        --- test orthonormality
c            ylmerr is the largest deviation from orthonormality
c            and will be used as a measure of the "machine zero".
         call tstylm(lmax,nomega,nlm,itch(ylm),itch(angwt),
     $               aitch(ptlm),itch(scr),ylmerr)
         write(iout,*) 'maximum ylm error',ylmerr
c
c        --- put down radial quadrature
         if(rtox) then
            rpts=ctheta
            rwt=rpts+nradial
            nrings=(nradial/(ncrule-1))+1
            jacob=rwt+nrings*(ncrule*(ncrule-1))
            npring=wpadti(jacob+nradial)
            top=iadtwp(npring+nradial)
            call cotes(alpha(ian(iatom)),nradial,ncrule,nrings,
     $                 itch(rpts),itch(rwt),itch(jacob),
     $                 aitch(npring),nrings,nr,nwtot)
         else
            rpts=ctheta
            rwt=rpts+mxrings*mxcotes
            npring=wpadti(rwt+mxrings*mxcotes*(mxrings*mxcotes-1))
            scr=iadtwp(npring+mxrings)
            top=scr+mxrings+1
            call newton(ian(iatom),nrings,itch(rpts),itch(rwt),
     $                  itch(scr),aitch(npring),mxrings,nr,nwtot)
         endif
         if(debug) then
            write(iout,*) 'nrings,nr',nrings,nr
            write(iout,*) 'n per ring',(aitch(npring+k),k=0,nrings-1)
            write(iout,*) 'rpts',(itch(rpts+k),k=0,nr-1)
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
         ngrid=nr*nomega
         gridwts=grid+3*ngrid
         ptrad=wpadti(gridwts+ngrid)
         top=iadtwp(ptrad+nr+1)
         call catenat(nr,ngrid,nomega,nomega,itch(rpts),itch(rwt),
     $               itch(angpts),itch(angwt),itch(grid),
     $               itch(gridwts),aitch(ptrad),c(1,iatom))
c
c        --- compute the atomic block analytically.
         va=gridwts
         vel=va+ngrid
         vnuc=vel+nnp*ngrid 
         nugrid=vnuc+nat
         scr=nugrid+3*ngrid
         call getscm(0,itch(scr),maxcor,'gofish',0)
c        the routine vatom expects the grid to be (3,ngrid)
c        at this point we have grid(ngrid,3)
c        turn it around.
         k=0
         do 41 pt=1,ngrid
            do 40 coord=1,3
               itch(nugrid+k)=itch(grid+(coord-1)*ngrid+pt-1)
               k=k+1
   40       continue
   41    continue
         call vatom(c,ex,itch(scr),cont,ptprim,noprim,nocont,ptcont,
     $           nat,nprim,iadtwp(maxcor),ntypes,nbtype,nnp,ncont,
     $           start,nbf,itch(zan),nocart,nobf,maxmom,mintyp,
     $           nx,ny,nz,minmom,' ',itch(vel),itch(vnuc),.false.,
     $           ' ',ngrid,itch(nugrid),d,iatom,itch(va))
         if(debug) then
            write(iout,*) 'atomic analytic'
            write(iout,*) (itch(va+i-1),i=1,ngrid)
         endif
c
c        --- decompose atomic potential into ylm components
ctmp     construct a xe-r2 potential
ctmp        write(iout,*) 'point,potne'
ctmp         do 311 i=1,ngrid
ctmp            ptx=grid+i-1
ctmp            pty=grid+ngrid+i-1
ctmp            ptx=grid+i-1
ctmp            ptz=grid+2*ngrid+i-1
ctmp            r2=itch(ptx)*itch(ptx)+itch(pty)*itch(pty)
ctmp     $        +itch(ptz)*itch(ptz)
ctmp            itch(va+i-1)=itch(ptx)*exp(-r2)
ctmp            write(iout,*) itch(ptx),itch(pty),itch(ptz)
ctmp            write(iout,*) itch(va+i-1)
ctmp  311    continue 
ctmp         call vmove(itch(poten),itch(va+nomega),ngrid-nomega)

         valm=vel
         scr=valm+nr*nlm
         top=scr+nr*nlm
         call ftoylm(itch(va),itch(valm),itch(ylm),
     $               nr,nomega,nlm,itch(angwt),itch(scr),ylmerr)
         if(debug) then
            write(iout,*) 'valm'
            do 210 i=1,nlm
               write(iout,*) 'lm',(itch(valm+(i-1)*nr+k),k=0,nr-1)
  210       continue
         endif
c
c
c        --- multiply va(r) by r to get something to spline fit.
         ioff=0
         falm=scr
         top=scr+nr*nlm
         do 212 i=1,nlm
            call vmul(itch(falm+ioff),itch(valm+ioff),itch(rpts),nr)
            if(debug) then
               write(iout,*) 'falm',(itch(falm+ioff+k),k=0,nr-1)
            endif
            ioff=ioff+nr
  212    continue
c
c        --- find voronoi wts for this atom
         vwts=top
         wts=vwts+ngrid
         rnuc=wts+ngrid
         amu=rnuc+nat*nat
         pwtx=amu+nat*nat
         rr=pwtx+nat
         top=rr+nat
         call voronoi(nat,c,itch(grid),itch(wts),ngrid,ngrid,
     $                iatom,itch(vwts),itch(rnuc),itch(amu),
     $                itch(pwtx),itch(rr),adjust,itch(radii))
         if(debug) then
            write(iout,*) 'vwts',(itch(vwts+k),k=0,ngrid-1)
         endif
c        --- put down the basis functions on the grid.
         phibf=wts
         grad=wts
         scrgrd=phibf+3*ngrid*nbf
         rsq=scrgrd+ngrid
         s=rsq+ngrid
         r=s+ngrid*mxcont
         xyzpow=r+ngrid*mxcont
         top=xyzpow+ngrid*3*max(bigl,2)
         ioff=0
c        assume there is enough left
         call bfgrd(c,ex,cont,ptprim,noprim,nocont,ptcont,
     $              nat,nprim,ntypes,nbtype,nnp,ncont,
     $              start,nbf,nocart,nobf,maxmom,minmom,mintyp,
     $              nx,ny,nz,itch(grid),ngrid,
     $              itch(xyzpow),itch(scrgrd),itch(rsq),itch(s),
     $              itch(r),itch(phibf),itch(grad),dograd,mxcont,
     $              ngrid,ngrid,ioff,maxl)
c
c        --- form density on the grid for this atom,
c            must worry about ndmat here someday
         dengrid=rsq
         dengrad=rsq
         top=dengrid+ngrid
         call rzero(itch(dengrid),ngrid)
c        must modify the one below to do d(1,dmat)
c        modify the density matrix to omit the diagonal blocks.
         call gridden(nbf,nnp,ngrid,ngrid,itch(d2),itch(phibf),
     $                itch(grad),itch(dengrid),itch(dengrad),
     $                itch(scrgrd),dmcut,.false.)
         if(debug) then
            write(iout,*) 'dengrid',(itch(dengrid+k),k=0,ngrid-1)
         endif
c
c        --- multiply by the voronoi weights to get the atomic
c            component of the density.
         call vmul(itch(dengrid),itch(dengrid),itch(vwts),ngrid)
         if(debug) then
            write(iout,*) 'weighted',(itch(dengrid+k),k=0,ngrid-1)
         endif
c
c        --- decompose density into ylm components
         rholm=top
         scr=rholm+nr*nlm
         top=scr+nr*nlm
         call ftoylm(itch(dengrid),itch(rholm),itch(ylm),
     $               nr,nomega,nlm,itch(angwt),itch(scr),ylmerr)
         if(debug) then
            write(iout,*) 'flm'
            do 111 i=1,nlm
               write(iout,*) 'lm',(itch(rholm+(i-1)*nr+k),k=0,nr-1)
  111       continue
         endif
c
c        --- solve the radial equations
         vlm=scr
         int0f=vlm+nr*nlm
         int0b=int0f+nlm
         scr=int0b+nlm
         tmp=scr+nr
         j=tmp+nr
         y=j+nr*(lmax+1)
         top=y+nr*(lmax+1)
         call rsolver(nrings,nr,aitch(npring),lmax,nlm,
     $                aitch(ptlm),itch(rpts),itch(rwt),nwtot,
     $                itch(rholm),itch(j),itch(y),itch(scr),
     $                itch(tmp),itch(vlm),itch(int0f),itch(int0b))
         if(debug) then
            write(iout,*) 'iatom,int0f',iatom,(itch(int0f+i-1),i=1,nlm)
            write(iout,*) 'iatom,int0b',iatom,(itch(int0b+i-1),i=1,nlm)
            ioff=0
            do 71 lm=1,nlm
               call vadd(itch(vlm+ioff),
     $                   itch(vlm+ioff),itch(falm+ioff),nr)
               write(iout,*) 'ulm diffs',lm
                  do 72 i=1,nr
                     diff=abs(itch(falm+ioff+i-1)-itch(vlm+ioff+i-1))
                     if(diff.gt.1.0d-08) then
c                       write(iout,*) 'pt,diff',i,diff
                     endif
   72             continue
               ioff=ioff+nr
   71       continue
         endif
c 
c        --- spline fit the radial solutions
         sc=scr 
         top=sc+(4*nr+1)*order
         call prespl(nr,itch(rpts),order,itch(sc))
         call splmat(nr,itch(rpts),order,itch(sc))
c
         coef=top
         break=coef+order*nr*nlm 
         top=break+nr+1
c        --- don't forget to combine the analytical and numerical radial
c            contributions
         ioff=0
         joff=0
         call rzero(itch(coef),order*nr*nlm)
         do 60 lm=1,nlm
            if(debug) then
               write(iout,*) 'vlm',(itch(vlm+ioff+i-1),i=1,nr)
            endif
ctmp            call vadd(itch(vlm+ioff),itch(vlm+ioff),itch(falm+ioff),nr)
            call vmove(itch(vlm+ioff),itch(falm+ioff),nr)
            if(debug) then
               write(iout,*) 'sum to fit',lm
               do 61 i=1,nr
                  write(iout,*) itch(rpts+i-1),itch(vlm+ioff+i-1)
   61          continue
            endif
            call splcof(nr,itch(vlm+ioff),order,nbreak,itch(break),
     $                  itch(coef+joff),itch(sc))
            write(iout,*) 'nbreak',nbreak
            ioff=ioff+nr
            joff=joff+nbreak*order
   60    continue
c
c        --- put the solution back onto the full grid
c            note that here we use the original grid.
         do 80 jatom=1,nat
            joff=0
            do 70 jgblk=1,ngb(jatom)
               ngj=gblksz(jgblk,jatom)
               fr=top
               ylm=fr+ngj*nlm
               ptlm=wpadti(ylm+ngj*nlm)
               ctheta=iadtwp(ptlm+(lmax+1)*(2*lmax+1))
               phi=ctheta+ngj
               cphi=phi+ngj
               sphi=cphi+ngj
               plm=sphi+ngj
               scr=plm+ngj*(lmax+1)
               ind=wpadti(scr+ngj)
c              top=ind+ngj
c
               ptv=joff+mxgrd*(jatom-1)
               call ftogrid(ngj,lmax,nlm,mxgrd,joff,
     $                      xyzgrid(1,1,jatom),c(1,iatom),itch(v+ptv),
     $                      itch(fr),itch(ylm),aitch(ptlm),itch(plm),
     $                      itch(ctheta),itch(phi),itch(cphi),
     $                      itch(sphi),itch(scr),nbreak,order,
     $                      itch(coef),itch(break),aitch(ind),
     $                      alpha(ian(iatom)))
               if(debug) then
                  write(iout,*) 'fitted v'
                  write(iout,*) (itch(v+ptv+i-1),i=1,ngj)
               endif
               joff=joff+ngj
   70       continue
   80    continue
c
c
  100 continue
c
c     --- compare the result with analytical one
      if(compare) then
         write(iout,*) 
     $      'compare:poisson-analytical>1.0d-06'
         mxerr=zero
         terr=zero
         ntest=0
         do 110 iatom=1,nat
            ioff=(iatom-1)*mxgrd
            pt=ioff
            do 109 igblk=1,ngb(iatom)
               ngi=gblksz(igblk,iatom)
               do 108 i=1,ngi
                  diff=itch(v+pt+i-1)-itch(poten+pt+i-1)
                  mxerr=max(mxerr,abs(diff))
                  terr=terr+diff
                  ntest=ntest+1
                   if(abs(diff).ge.1.0d-03) then
                      write(iout,*) 'point,diff',
     $                        i,diff,itch(v+pt+i-1),itch(poten+pt+i-1)
                      write(iout,*) xyzgrid(i,1,iatom),
     $                            xyzgrid(i,2,iatom),xyzgrid(i,3,iatom)
                   endif
  108          continue
               pt=pt+ngi
  109       continue
  110    continue
         write(iout,*) 'maximum error',mxerr
         write(iout,*) 'average error',terr/ntest
      endif
c
c     --- at this stage we have the coulomb kernel in the array v.
c         will want to put it on the grid now to get jmat.
      return
      end
