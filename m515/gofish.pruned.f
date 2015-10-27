*deck @(#)gofish.pruned.f	5.1  11/28/95
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
c***begin prologue     gofish.pruned.f
c***date written       940304      (yymmdd)  
c***revision date      11/28/95
c
c***keywords           poisson, coulomb, potential, density
c***author             martin, richard(lanl) 
c***source             @(#)gofish.pruned.f	5.1 11/28/95
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
c***end prologue       gofish.pruned.f
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
c THIS STUFF WON'T WORK WITH PRUNED GRIDS, need to get ylms for lots of 
C  different lebedev grids instead.  Use grid file I/O and not on-the-fly
c  calculations.
c$$$c     --- generate the angular grid for decomposition.
c$$$c            lebedev angular points and weights on the 
c$$$c            unit sphere.
c$$$      lebord=2*vlmax+1
c$$$      vnomega=angsiz(lebord)
c$$$c     --- allocate some core. these first few will contain the radial
c$$$c         points, the potential, and the spline coefficients.
c$$$      angpts=1
c$$$      angwt=angpts+3*vnomega
c$$$      top0=angwt+3*vnomega
c$$$      if (top0 .gt. left) then
c$$$         write(iout,*)'not enough core for lebedev !'
c$$$         write(iout,*)'have ',left,' need ',top0
c$$$         call lnkerr('m511: gofish')
c$$$      endif
c$$$      if(timeit) then
c$$$         call timing(dum1,dum2,dum3)
c$$$      endif
c$$$      call sphere(vnomega,itch(angpts),itch(angwt),lebord,vnomega)
c$$$      if(debug) then
c$$$         write(iout,*) 'angtps',(itch(angpts+k),k=0,3*vnomega-1)
c$$$         write(iout,*) 'angwt',(itch(angwt+k),k=0,vnomega-1)
c$$$      endif
c$$$c
c$$$c     --- generate ylm on the lebedev grid
c$$$      ylm=top0
c$$$      ptlm=wpadti(ylm+vnomega*nlm)
c$$$      ctheta=iadtwp(ptlm+(vlmax+1)*(2*vlmax+1))
c$$$      phi=ctheta+vnomega
c$$$      cphi=phi+vnomega
c$$$      sphi=cphi+vnomega
c$$$      plm=sphi+vnomega
c$$$      scr=plm+vnomega*(vlmax+1)
c$$$      top0=scr+vnomega
c$$$      if (top0 .gt. left) then
c$$$         write(iout,*)'not enough core for vylm !'
c$$$         write(iout,*)'have ',left,' need ',top0
c$$$         call lnkerr('m511: gofish')
c$$$      endif
c$$$      ioff=0
c$$$      call vylm(vlmax,vnomega,ioff,itch(angpts),test,vnomega,nlm,
c$$$     $          itch(ylm),aitch(ptlm),itch(ctheta),itch(phi),
c$$$     $          itch(cphi),itch(sphi),itch(plm),itch(scr))
c$$$      if(debug) then
c$$$         do 113 i=1,nlm
c$$$            write(iout,*) 'ylm;l,m',i,
c$$$     $         (itch(ylm+(i-1)*vnomega+k-1),k=1,vnomega)
c$$$  113    continue
c$$$      endif
c$$$c
c$$$c     --- test orthonormality
c$$$c         ylmerr is the largest deviation from orthonormality
c$$$c         and will be used as a measure of the "machine zero".
c$$$      call tstylm(vlmax,vnomega,nlm,itch(ylm),itch(angwt),
c$$$     $            aitch(ptlm),itch(scr),ylmerr)
c$$$c      write(iout,*) 'maximum ylm error',ylmerr
c$$$      if(timeit) then
c$$$         call timing(dum4,dum5,dum6)
c$$$         timpre=timpre+dum4-dum1
c$$$      endif
c
c
c
c     --- generate the potential originating from this charge 
c         distribution.
c$$$      if(dogrdio) then
c$$$         if (.not.called) then
c$$$            call iosys('create real "poisson voronoi weights" on '
c$$$     $           //grdfil(1:3),pmxgrd*nat,0,0,' ')
c$$$         endif
c$$$c         call iosys('rewind "poisson voronoi weights" on '//grdfil(1:3),
c$$$c     $              0,0,0,' ')
c$$$      endif
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
c        --- get this atom's radial grid, and info about its angular grids
c            from the grid file
c        --- get this atom's combined radial/angular grid from grid file
c        --- get this atom's voronoi weights from grid file
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
C NEW CODE: 
c    zero all rholm channels
C    for each block of radial shells on this atom with like angular grids,
c    read in the appropriate set of ylms, fit each radial shell in this
c    block to ylms up to the appropriate order.
c    Some rholms will wind up being zero at spots because of the differing
c    orders.

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
         call rsolver(nrings,nr,aitch(npring),vlmax,nlm,
     $                aitch(ptlm),rpts(1,iatom),itch(rwt),nwtot,
     $                itch(rholm),itch(j),itch(y),itch(scr),
     $                itch(tmp),vlm(1,1,iatom),itch(int0f),itch(int0b))
c$$$         write(iout,*)' vlms'
c$$$         do 666 i=1,nlm
c$$$            write(iout,*)' channel ',i
c$$$            write(iout,*)(vlm(ig,i,iatom),ig=1,nr)
c$$$ 666     continue 
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
         ioff=0
         do 60 lm=1,nlm
            if(debug) then
               write(iout,*) 'vlm',(vlm(i,lm,iatom),i=1,nr)
            endif
            if(debug) then
              write(iout,*) 'sum to fit',lm
              do 61 i=1,nr
                 write(iout,*) rpts(i,iatom),vlm(i,lm,iatom)
   61         continue
            endif
            call spline3(rpts(1,iatom),vlm(1,lm,iatom),nr,1.0d+31,
     $                   1.0d+31,y2(1,lm,iatom),itch(scr))
            ioff=ioff+nr
   60    continue
         if (dumpfit) then
            do 6060 l=0,vlmax
               do 6061 m=-l,l
                  write(4,*)'atom, l,m',iatom,l,m
                  write(4,*)'sum to fit'
                  do 6062 i=1,nr
                     write(4,*) rpts(i,iatom),vlm(i,
     $                    aitch(ptlm+l+(m+vlmax)*(vlmax+1)),iatom)
 6062             continue 
 6061          continue 
 6060       continue 
         endif
         if(timeit) then
            call timing(dum4,dum5,dum6)
            timrsolv=timrsolv+dum4-dum1
         endif
         next=nxtval(nproc)+1
   80 continue
      next=nxtval(-nproc)
c spread the spline data across the nodes --- do this with a global sum since
c all the stuff we didn't calculate is zeroed out anyway.
      if (ispar) then
         call dgop(601,rpts,vradial*nat,'+')
         call dgop(602,vlm,vradial*nlm*nat,'+')
         call dgop(603,y2,vradial*nlm*nat,'+')
      endif
c
c     --- we now have the potential decomposed into atomic contributions
c         and fit to a cubic spline: rpts,vlm,y2.
c
c         --- put the solution back onto the full grid
c             note that here we use the original grid.
      top0=1
      if(dogrdio) then
         grdsiz=1
         top0=iadtwp(grdsiz+nat)
         call iosys('read integer "external atomic grid size" from '
     $        //grdfil(1:3),nat,aitch(grdsiz),0,' ') 
      endif
      next=nxtval(nproc)+1
      do 100 jatom=1,nat
         if (jatom .ne. next) goto 100
         if (grdtyp(jatom).eq.'nobasis') goto 100
c
c        --- generate grid, calc number of grid blocks
         call timing(dum1,dum2,dum3)
         if(dogrdio) then
            jngrid=aitch(grdsiz+jatom-1)
            call iosys('read real "external grid" from '
     $           //grdfil(1:3),mxgrd*3,xyzgrid,mxgrd*3*(jatom-1),' ')
            call iosys('read real "external grid weights" from '
     $                 //grdfil(1:3),mxgrd,grdwts,mxgrd*(jatom-1),' ')
         else
            call mkatmg(c,ian,xyzgrid,grdwts,rmax,lmax,nomega,
     $                  nradial,
     $                  nat,jngrid,mxgrd,vwts,rnuc,amu,pwtx,rr,adjust,
     $                  radii,akl,grdtyp(jatom),
     $                  radshls,ptrad,.false.,vwts,jatom)
         endif
c         if(timeit) then
c            call timing(dum4,dum5,dum6)
c            timatg=timatg+dum4-dum1
c         endif
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
            r=fr+ngj*nlm
            ylm=r+ngj
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
               call lnkerr('m511: gofish')
            endif
c
            call rzero(itch(v),ngj)
            if(timeit) then
               call timing(dum1,dum2,dum3)
            endif
            do 75 iatom=1,nat
               if (grdtyp(iatom).eq.'nogrid') goto 75
               if(method1) then
                  call ftogrid(ngj,vlmax,nlm,mxgrd,joff,
     $                         xyzgrid,c(1,iatom),itch(v),
     $                         itch(fr),itch(ylm),aitch(ptlm),
     $                         itch(plm),itch(ctheta),itch(phi),
     $                         itch(cphi),
     $                         itch(sphi),itch(scr),nr,rpts(1,iatom),
     $                         vlm(1,1,iatom),y2(1,1,iatom),aitch(ind))
               else
                  call ftogrid2(ngj,minesz,vlmax,nlm,mxgrd,joff,
     $                         xyzgrid,c(1,iatom),itch(v),
     $                         itch(fr),itch(r),itch(ylm),aitch(ptlm),
     $                         itch(plm),itch(ctheta),itch(phi),
     $                         itch(cphi),
     $                         itch(sphi),itch(scr),nr,rpts(1,iatom),
     $                         vlm(1,1,iatom),y2(1,1,iatom),aitch(ind))
               endif
   75       continue
c            if(timeit) then
c               call timing(dum4,dum5,dum6)
c               tim2grd=tim2grd+dum4-dum1
c            endif
            call timing(dum4,dum5,dum6)
            tim2grd=tim2grd+dum4-dum1
            call timing(dum1,dum2,dum3)
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
     $                 mxgbsiz,ngj,joff,maxl,bfcut)
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
     $                   jmat,kmcut)
            joff=joff+ngj
c            if(timeit) then
c               call timing(dum4,dum5,dum6)
c               timfmk=timfmk+dum4-dum1
c            endif
            call timing(dum4,dum5,dum6)
            timfmk=timfmk+dum4-dum1
   90    continue
         next=nxtval(nproc)+1
  100 continue
      next=nxtval(-nproc)
c
c     --- as a temp fix to the three center stuff, now zero
c         out any piece of the jmatrix which corresponds to two functions
c         on the same site. 
      if(rmdiag) then
         d2=1
         top0=d2+nnp
         if (top0 .gt. left) then
            write(iout,*)'not enough core for d2!'
            write(iout,*)'have ',left,' need ',top0
            call lnkerr('m511: gofish')
         endif
         call vmove(itch(d2),jmat,nnp)
         do 49 iatom=1,nat
            jatom=iatom
            do 47 itype=1,nbtype
               if (noprim(iatom,itype).le.0) go to 47
               if (iatom.ne.jatom) then
                  jtypmx=nbtype
               else
                  jtypmx=itype
               end if
               do 46 jtype=1,jtypmx
                  if (noprim(jatom,jtype).le.0) go to 46
c
                  nconti=nocont(iatom,itype)
                  ncontj=nocont(jatom,jtype)
                  lenblk=nocart(itype)*nocart(jtype)
                  call rzero(itch(top0),nconti*ncontj*lenblk)
                  call put1el(itch(d2),itch(top0),start,iatom,jatom,
     $                        itype,jtype,nconti,ncontj,nnp,lenblk,nat,
     $                        nbtype,nobf)
   46          continue
   47       continue
   49    continue
   50    continue
         call vmove (jmat,itch(d2),nnp)
         if(debug) then
            write(iout,*) 'nnp',nnp
            write(iout,*) 'jmat',(jmat(i,1),i=1,nnp)
         endif
      endif
c   41 continue
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
