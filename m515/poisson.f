*deck @(#)poisson.f	5.1  11/28/95
      subroutine poisson(values,d,dlast,nbf,nnp,jmat,t1,ncoul,
     $                   ndmat,nat,mxgrd,pmxgrd,
     $                   dmcut,dencut,bfcut,kmcut,
     $                   mxgbsiz,valuesi,ian,
     $                   coord,ex,cont,ptprim,noprim,nocont,ptcont,
     $                   mxcont,nprim,ntypes,nbtype,ncont,
     $                   start,nocart,nobf,maxmom,minmom,mintyp,
     $                   nx,ny,nz,lenxyz,xyzgrid,grdwts,charge,bigl,
     $                   ops,left,vwts,rnuc,amu,pwtx,rr,radii,akl,ptrad,
     $                   rmax,lmax,nomega,nradial,grdtyp,adjust,minesz,
     $                   vlmax,vradial,vncrule,mxiter,ntotal,nonzer,npf,
     $                   nnprim,jprim,nnshl,dijmax,qint,qtest,cutexp,
     $                   rhotest,prtoao,pstart,dpr,method,grdfil)
c***begin prologue     poisson.f
c***date written       940304      (yymmdd)  
c***revision date      11/28/95
c
c***keywords           poisson, coulomb, potential, density
c***author             martin, richard(lanl) 
c***source             @(#)poisson.f	5.1   11/28/95
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
c***end prologue       %%
      implicit none
c     --- input variables ---
      integer nbf,nnp,ncoul,ndmat,nat,mxgrd,pmxgrd,mxgbsiz
      integer nprim,ntypes,nbtype,ncont,mxcont
      integer bigl
      integer left
      integer rmax,lmax,nomega,nradial,minesz
      integer nlm,vlmax,vradial,vncrule
      integer mxiter,npf,nnprim,nnshl,lenxyz
      integer ntotal,nonzer
      real*8 dencut,dmcut,bfcut,kmcut
      real*8 qtest,cutexp
      logical adjust,rhotest
      character*(*) method
      character*4 grdfil
c     --- input arrays (unmodified) ---
      character*(*) ops
      character*(*) grdtyp(nat)
      integer ian(nat)
      integer ptprim(nat,ntypes),noprim(nat,ntypes),nocont(nat,ntypes)
      integer ptcont(nat,ntypes),start(nat,ntypes),nocart(0:*)
      integer nobf(ntypes),maxmom(ntypes),minmom(ntypes),mintyp(ntypes)
      integer nx(lenxyz),ny(lenxyz),nz(lenxyz)
      integer prtoao(npf,nbf),pstart(nat,nbtype)
      real*8 d(nnp,ndmat),dlast(nnp,ndmat)
      real*8 coord(3,nat),ex(nprim),cont(ncont)
      real*8 xyzgrid(mxgrd,3),grdwts(mxgrd,nat)
c     --- input arrays (scratch) ---
      real*8 vwts(mxgrd),rnuc(nat,nat),amu(nat,nat),pwtx(nat)
      real*8 rr(nat),radii(nat),akl(nat,nat)
      real*8 values(left),t1(nnp,ncoul)
      real*8 jprim(nnprim),dpr(nnprim,ndmat),qint(nnshl),dijmax(nnshl)
      integer valuesi(*)
      integer ptrad(rmax)
c     --- output arrays ---
      real*8 jmat(nnp,ncoul)
      real*8 charge(nat,ndmat)
c     --- output variables ---
c     --- scratch arrays ---
c     --- local variables ---
      integer inp,iout
      integer iadtwp,wpadti
      integer dtot,i
      integer imaxl,itch,aitch,pleft,rpts,vlm,y2
      logical do0,do2,do3
      logical zdiag,zoffdiag,rmdiag
      logical prnt,debug
      logical timeit
      real*8 dum1,dum2,dum3,dum4,dum5,dum6
      real*8 timej(3)
      data prnt/.true./
      parameter (debug=.false.,timeit=.false.)
      integer stderr,ierr
c
      common/io/inp,iout
c      
      integer mynodeid, nodeid
      logical ispar
      common /tcgmesa/ ispar

      ierr=stderr()
      mynodeid=nodeid()
c
c     --- evaluate the j-matrix using a combination of analytical integrals
c         and the poisson solver.
c
c     --- some core necessary for gofish.
      dtot=1
      imaxl=wpadti(dtot+nnp)
      itch=iadtwp(imaxl+nat)
      aitch=wpadti(itch)
      call rzero(values(dtot),nnp)
      pleft=left-nnp-nat
c
      nlm=(vlmax+1)*(vlmax+1)
      rpts=itch
      vlm=rpts+(vradial)*nat
      y2=vlm+(vradial)*nlm*nat
      itch=y2+(vradial)*nlm*nat
      aitch=wpadti(itch)
      pleft=pleft-(vradial)*nat*(1+nlm+nlm)
c
c     branch depending on the number of analytic ntegrals to be 
c     processed.
      do0=.false.
      do3=.false.
      do2=.false.
      if(method.eq.'do0') then
         do0=.true.
      else if(method.eq.'do3') then
         do3=.true.
      else if(method.eq.'do2') then
         do2=.true.
      else
         call lnkerr('unrecognized method in poisson')
      endif
      write(ierr,*) 'do0,do2,do3',do0,do2,do3
      if(do0) then
         write(ierr,*) 'Poisson: method do0'
c        --- do all the integrals via poisson.
         call rzero(t1,nnp)
         call rzero(jmat,nnp)
         ntotal=nnshl
c        --- pass the difference density to gofish.
         call vsub(values(dtot),d(1,1),dlast(1,1),nnp)
         zdiag=.false.
         zoffdiag=.false.
         rmdiag=.false.
         call timing(dum1,dum2,dum3)
         call gofish(values(itch),values(dtot),nbf,nnp,t1,ncoul,
     $               ndmat,nat,mxgrd,pmxgrd,dmcut,dencut,bfcut,
     $               kmcut,mxgbsiz,valuesi(aitch),ian,coord,ex,cont,
     $               ptprim,noprim,
     $               nocont,ptcont,mxcont,nprim,ntypes,nbtype,ncont,
     $               start,nocart,nobf,maxmom,minmom,mintyp,nx,ny,nz,
     $               xyzgrid,grdwts,charge,valuesi(imaxl),bigl,ops,
     $               pleft,vwts,rnuc,amu,pwtx,rr,radii,akl,ptrad,rmax,
     $               lmax,nomega,nradial,grdtyp,adjust,minesz,vlmax,
     $               vradial,vncrule,nlm,values(rpts),values(vlm),
     $               values(y2),zdiag,zoffdiag,rmdiag,grdfil)
         call vneg(jmat,t1,nnp)
         call timing(dum4,dum5,dum6)
         write(4,*)'    node ',mynodeid,' do0 time for gofish',
     $        dum4-dum1,dum5-dum2,dum6-dum1
      else if(do3) then
         call rzero(jmat,nnp*ncoul)
c           --- do only the 4-center integrals via poisson.
         call timing(dum1,dum2,dum3)
c           -- form the difference density matrix.
c              use t1 to hold the difference. note that it is given room
c              equivalent to max(nbf**2,nnp*ncoul) in the calling routine.
c              the density difference in the primitive basis comes back in dpr.
         do 132 i=1,ndmat
            call vsub(t1,d(1,i),dlast(1,i),nnp)
            call tr1dm(values(1),values(1+npf*npf),t1,
     $           prtoao,dpr(1,i),nbf,nnp,npf,nnprim)
            if(debug) then
               write(iout,*) 'difference primitive density matrices'
               call print(dpr(1,i),nnprim,npf,iout)
            endif
c
c              --- now form an array which contains the largest difference
c                  density in the contracted basis for each shell block
c                  combination. this will be used in conjunction with the
c                  largest estimated integral in the block to screen zeroes.
            call bigdij(noprim,nbtype,nocont,ncont,nocart,nat,
     $           start,nbf,nnp,t1,nnshl,dijmax,values,left,
     $           ndmat)
 132     continue
c           ---note that this call to directj requests only (ii,kl) and
c              (ij,kk) integrals be done. i.e. the 3-center ones.
         call directj(jmat,jprim,ptprim,noprim,nbtype,ex,coord,
     $        nx,ny,nz,lenxyz,nocart,mintyp,maxmom,
     $        values,left,nat,npf,nnprim,nprim,ops,cutexp,
     $        rhotest,values,.false.,ndmat,
     $        pstart,prtoao,dpr,
     $        nbf,nnp,ntotal,nonzer,nnshl,qint,
     $        dijmax,qtest,'closed',.false.,values,values,
     $        .true.,.false.)
c           --- add the 3-center contribution difference j-matrix 
c               to the last one ---
         call timing(dum4,dum5,dum6)
         timej(1)=timej(1)+dum4-dum1
         timej(2)=timej(2)+dum5-dum2
         write(4,*) '    node',mynodeid,
     $        ' do3time for 3-site j-matrix',dum4-dum1,dum5-dum2,
     $        dum6-dum3
c     
c
         call timing(dum1,dum2,dum3)
c           put the 4-center piece in t1.
         call rzero(t1,nnp)
c           --- pass the difference density to gofish.
         call vsub(values(dtot),d(1,1),dlast(1,1),nnp)
         zdiag=.true.
         zoffdiag=.false.
         rmdiag=.true.
         call gofish(values(itch),values(dtot),nbf,nnp,t1,ncoul,
     $        ndmat,nat,mxgrd,pmxgrd,dmcut,dencut,bfcut,kmcut,
     $        mxgbsiz,valuesi(aitch),ian,coord,ex,cont,ptprim,noprim,
     $        nocont,ptcont,mxcont,nprim,ntypes,nbtype,ncont,
     $        start,nocart,nobf,maxmom,minmom,mintyp,nx,ny,nz,
     $        xyzgrid,grdwts,charge,valuesi(imaxl),bigl,ops,
     $        pleft,vwts,rnuc,amu,pwtx,rr,radii,akl,ptrad,rmax,lmax,
     $        nomega,nradial,grdtyp,adjust,minesz,vlmax,vradial,
     $        vncrule,nlm,values(rpts),values(vlm),values(y2),
     $        zdiag,zoffdiag,rmdiag,grdfil)
c     add it to the 3-center piece.
         call vneg(t1,t1,nnp)
         call vadd(jmat,jmat,t1,nnp)
         call timing(dum4,dum5,dum6)
         timej(1)=timej(1)+dum4-dum1
         timej(2)=timej(2)+dum5-dum2
         write(4,*)'    node ',mynodeid,' time for gofish',
     $        dum4-dum1,dum5-dum2,dum6-dum3
         if(timeit) then
            write(iout,*) 'time for poisson',dum4-dum1
            write(iout,*) 'total time for j',timej(1)
         endif
      else if(do2) then
         call timing(dum1,dum2,dum3)
c           -- form the difference density matrix.
c              use t1 to hold the difference. note that it is given room
c              equivalent to max(nbf**2,nnp*ncoul) in the calling routine.
c              the density difference in the primitive basis comes back in dpr.
         do 232 i=1,ndmat
            call vsub(t1,d(1,i),dlast(1,i),nnp)
            call tr1dm(values(1),values(1+npf*npf),t1,
     $           prtoao,dpr(1,i),nbf,nnp,npf,nnprim)
            if(debug) then
               write(iout,*) 'difference primitive density matrices'
               call print(dpr(1,i),nnprim,npf,iout)
            endif
c
c              --- now form an array which contains the largest difference
c                  density in the contracted basis for each shell block
c                  combination. this will be used in conjunction with the
c                  largest estimated integral in the block to screen zeroes.
            call bigdij(noprim,nbtype,nocont,ncont,nocart,nat,
     $           start,nbf,nnp,t1,nnshl,dijmax,values,left,
     $           ndmat)
 232     continue
c           ---note that this call to directj requests only (ii,kk)
c              integrals be done, i.e. the 2-center ones.
         call directj(jmat,jprim,ptprim,noprim,nbtype,ex,coord,
     $        nx,ny,nz,lenxyz,nocart,mintyp,maxmom,
     $        values,left,nat,npf,nnprim,nprim,ops,cutexp,
     $        rhotest,values,.false.,ndmat,
     $        pstart,prtoao,dpr,
     $        nbf,nnp,ntotal,nonzer,nnshl,qint,
     $        dijmax,qtest,'closed',.false.,values,values,
     $        .true.,.true.)
         call timing(dum4,dum5,dum6)
         timej(1)=timej(1)+dum4-dum1
         timej(2)=timej(2)+dum5-dum2
         write(4,*)'    node ',mynodeid,' do2 time for jmatrix',
     $        dum4-dum1,dum5-dum2,dum6-dum3
         if(timeit) then
            write(iout,*) 'time for j-matrix',
     $           dum4-dum1,dum5-dum2,dum6-dum3
         endif
c     
c           --- call it the first time with just off-diagonals
         call timing(dum1,dum2,dum3)
c           put the 4-center piece in t1.
         call rzero(t1,nnp)
c           --- pass the difference density to gofish.
         call vsub(values(dtot),d(1,1),dlast(1,1),nnp)
c           zero the diagonal and off-diagonal blocks, but don't
c           remove the diagonals.
         zdiag=.true.
         zoffdiag=.false.
         rmdiag=.false.
         call gofish(values(itch),values(dtot),nbf,nnp,t1,ncoul,
     $        ndmat,nat,mxgrd,pmxgrd,dmcut,dencut,bfcut,kmcut,
     $        mxgbsiz,valuesi(aitch),ian,coord,ex,cont,ptprim,noprim,
     $        nocont,ptcont,mxcont,nprim,ntypes,nbtype,ncont,
     $        start,nocart,nobf,maxmom,minmom,mintyp,nx,ny,nz,
     $        xyzgrid,grdwts,charge,valuesi(imaxl),bigl,ops,
     $        pleft,vwts,rnuc,amu,pwtx,rr,radii,akl,ptrad,rmax,lmax,
     $        nomega,nradial,grdtyp,adjust,minesz,vlmax,vradial,
     $        vncrule,nlm,values(rpts),values(vlm),values(y2),
     $        zdiag,zoffdiag,rmdiag,grdfil)
         call vneg(t1,t1,nnp)
c           add it to the 2-center piece.
         call vadd(jmat,jmat,t1,nnp)
         call timing(dum4,dum5,dum6)
         timej(1)=timej(1)+dum4-dum1
         timej(2)=timej(2)+dum5-dum2
         write(4,*)'    node ',mynodeid,' do2 time for o-d gofish',
     $        dum4-dum1,dum5-dum2,dum6-dum3
         if(timeit) then
            write(iout,*) 'time for poisson',dum4-dum1
            write(iout,*) 'total time for j',timej(1)
         endif
c
c           --- now do the diagonal pieces.
         call timing(dum1,dum2,dum3)
c           put the 4-center piece in t1.
         call rzero(t1,nnp)
c           --- pass the difference density to gofish.
         call vsub(values(dtot),d(1,1),dlast(1,1),nnp)
         zdiag=.false.
         zoffdiag=.true.
         rmdiag=.true.
         call gofish(values(itch),values(dtot),nbf,nnp,t1,ncoul,
     $        ndmat,nat,mxgrd,pmxgrd,dmcut,dencut,bfcut,kmcut,
     $        mxgbsiz,valuesi(aitch),ian,coord,ex,cont,ptprim,noprim,
     $        nocont,ptcont,mxcont,nprim,ntypes,nbtype,ncont,
     $        start,nocart,nobf,maxmom,minmom,mintyp,nx,ny,nz,
     $        xyzgrid,grdwts,charge,valuesi(imaxl),bigl,ops,
     $        pleft,vwts,rnuc,amu,pwtx,rr,radii,akl,ptrad,rmax,lmax,
     $        nomega,nradial,grdtyp,adjust,minesz,vlmax,vradial,
     $        vncrule,nlm,values(rpts),values(vlm),values(y2),
     $        zdiag,zoffdiag,rmdiag,grdfil)
         call vneg(t1,t1,nnp)
c           add it to the 3-center piece.
         call vadd(jmat,jmat,t1,nnp)
         call timing(dum4,dum5,dum6)
         timej(1)=timej(1)+dum4-dum1
         timej(2)=timej(2)+dum5-dum2
         write(4,*)'    node ',mynodeid,' do2 time for d gofish',
     $        dum4-dum1,dum5-dum2,dum6-dum3
         if(timeit) then
            write(iout,*) 'time for poisson2',dum4-dum1
            write(iout,*) 'total time for j',timej(1)
         endif
      endif
c
c
c now put all the j matrices together if we're parallel
c
      if (ispar) call dgop(615,jmat,nnp*ncoul,'+')
      
      return
      end
