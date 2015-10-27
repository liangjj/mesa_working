*deck %W%  %G%
      subroutine pm725(z,a)
c***begin prologue     pm725.f
c***date written       940513   (yymmdd)
c***revision date      4/17/95
c***keywords           m725, link 721, dft, integrals, derivatives, gradient
c***author             martin,richard and russo,thomas(lanl)
c***source             %W%   %G%
c***purpose            computes derivatives of the exchange-correlation
c                      energy
c***description
c***references
c
c***routines called
c***end prologue       pm725.f
      implicit none
c     --- input variables -----
c     --- input arrays (unmodified) ---
c     --- input arrays (scratch) ---
      integer a(*)
      real*8 z(*)
c     --- output arrays ---
c     --- output variables ---
c     --- scratch arrays ---
c     --- local variables ---
      logical prnt
      logical logkey
      logical slater,becke,lyp,vwn,cnull
      logical adjust
      logical prteexch,calce
      integer maxatm
      integer minesz
      parameter (maxatm=2000)
      character*16 grdtyp(maxatm)
      character*4096 ops
      character*8 prtflg
      character*8 calc,corrf,exchf
      integer inp,iout
      integer nat,nbf,nnp,ndmat,nderiv,nd1e,nd2e
      integer f,d,d2e,grad,dkay,zan,c,top
      integer nprim,ncont,ntypes,nbtype,lenxyz,mxcont
      integer maxl,ptprim,noprim,ptcont,nocont,start
      integer nocart,nobf,minmom,maxmom,mintyp,nx,ny,nz
      integer cont,ex
      integer intkey,wpadti,iadtwp
      integer multip,nae,nbe,nexch,mxgrd,rmax,lmax
      integer mxgblk,nomega,ian,ptrad,radshls
      integer grid,wts,vwts,rnuc,amu,pwtx,rr,radii,nradial,gradwts,akl
      integer canget,values,valuesi,mxgbsiz
      integer charge,bigl 
      integer angsiz
      integer left,dftabs,kleft,dftcore,dfttmp
      real*8 dmcut,defcut,dencut,toosmall
      real*8 exc,eexch,exc511
      real*8 fpkey
      real*8 two
c
      parameter (two=2.0d+00)
      parameter (defcut=1.0d-16,toosmall=1.0d-50)
      integer mynodeid,nodeid,mdtob,mitob
      include 'msgtypesf.h'
      logical ispar
      common /tcgmesa/ ispar
c
      common /io/     inp,iout
c
      data prnt/.true./
      save prnt
      calce=.true.
c
c
      mynodeid=nodeid()
c
 1000 format(1x,'m725: exchange-correlation derivatives')
 1010 format (5x,'the exchange-correlation contribution to the scf ',
     $        'gradients:')
 1020 format (5x,'the exchange-correlation contribution to the scf ',
     $       'force constants:')
 1030 format (5x,'the two-electron contribution to the gradients')
 1040 format (5x,'the two-electron contribution to the force'
     $        //' constants')
 1050 format (5x,'the total scf first derivatives:')
 1060 format (5x,'the force-constants less cphf contributions')
 1070 format(5x,'maximum block size',17x,i9)
 1080 format(8x,'grid size; atom',i3,':',16x,i6,2x,i3,' blocks')
c
      ispar=.true.
      if (mynodeid.eq.0) then
c     --- collect the options string ---
         call iosys('read character options from rwf',-1,0,0,ops)
         ispar=logkey(ops,'parallel',.false.,' ')
c
         call iosys('read character "print flag" from rwf',-1,0,0,
     $        prtflg)
         if(prtflg.eq.'minimum') prnt=.false.
c
         if(prnt) then
            write(iout,1000)
         endif
      endif
c
c     --- find out how much core is available ---
      call getscm(0,z,canget,'pm725',0)
c
c     --- set up some parameters depending on multip -----
      if (mynodeid.eq.0) then
         call iosys('read integer "spin multiplicity" from rwf',
     $        1,multip,0,' ')
         call iosys(
     $        'read integer "number of alpha electrons" from rwf',
     $        1,nae,0,' ')
         call iosys(
     $        'read integer "number of beta electrons" from rwf',
     $        1,nbe,0,' ')
      endif
      if (ispar) then
         call brdcst(1+MSGINT,multip,mitob(1),0)
         call brdcst(2+MSGINT,nae,mitob(1),0)
         call brdcst(3+MSGINT,nbe,mitob(1),0)
      endif
      if(multip.eq.1) then
         calc='closed'
         nexch=1
      else
         calc='open'
         nexch=2
      endif
c
c     --- pick a functional set (exchange,correlation)
      if (mynodeid.eq.0) then
         slater=logkey(ops,'scf=exchf=slater',.false.,' ')
         becke=logkey(ops,'scf=exchf=becke',.false.,' ')
         lyp=logkey(ops,'scf=corrf=lyp',.false.,' ')
         vwn=logkey(ops,'scf=corrf=vwn',.false.,' ')
         cnull=logkey(ops,'scf=corrf=null',.false.,' ')
         if (becke .and. slater)
     $        call plnkerr('m511: two correlation functionals chosen',1)
         if (lyp .and. vwn)
     $        call plnkerr('m511: two exchange functionals chosen',2)
c     --- have to pick at least one, make slater-vwn default
         if (.not.(becke .or. slater)) slater=.true.
         if (.not.(vwn.or.lyp.or.cnull)) vwn=.true.
      endif
      if (ispar) then
         call brdcst(4,slater,4,0)
         call brdcst(5,vwn,4,0)
         call brdcst(6,becke,4,0)
         call brdcst(7,lyp,4,0)
         call brdcst(8,cnull,4,0)
         if(slater) then
            exchf='slater'
         else if(becke) then
            exchf='becke '
         endif
         if(vwn) then
            corrf='vwn'
         else if(lyp) then
            corrf='lyp'
         else
            corrf='null'
         endif
      endif
c
c     --- quadrature options ---
      if (mynodeid.eq.0)
     $     call iosys('read character "atomic grid name" from rwf',
     $           -1,0,0,grdtyp)
      if (ispar) call brdcst(9,grdtyp,16*maxatm,0)
c     --- default mxgrd to the larges standard grid size now in effect.
c     this would be nr=50,nang=302
      if (mynodeid.eq.0) then
         mxgrd=intkey(ops,'scf=mxgrid',15100,' ')
         rmax=intkey(ops,'scf=radgrid',51,' ')
         lmax=intkey(ops,'scf=lebord',23,' ')
         dmcut=fpkey(ops,'scf=denmat-cutoff',defcut,' ')
         dencut=fpkey(ops,'scf=density-cutoff',toosmall,' ')
         minesz=intkey(ops,'scf=minesz',100,' ')
         mxgblk=intkey(ops,'scf=maxgblk',5,' ')
         adjust=logkey(ops,'scf=adjustcell',.false.,' ')
      endif
      if (ispar) then
         call brdcst(10+MSGINT,mxgrd,mitob(1),0)
         call brdcst(11+MSGINT,rmax,mitob(1),0)
         call brdcst(12+MSGINT,lmax,mitob(1),0)
         call brdcst(13+MSGDBL,dmcut,mdtob(1),0)
         call brdcst(14+MSGDBL,dencut,mdtob(1),0)
         call brdcst(15+MSGINT,minesz,mitob(1),0)
         call brdcst(16+MSGINT,mxgblk,mitob(1),0)
         call brdcst(17,adjust,4,0)
      endif
c
c     --- get the lengths of arrays needed for core allocation ---
      if (mynodeid .eq. 0) then
         call iosys('read integer "number of atoms" from rwf',1,nat,0,
     $        ' ')
         call iosys('read integer "number of basis functions" from rwf',
     $            1,nbf,0,' ')
      endif
      if(ispar) then
         call brdcst(18+MSGINT,nat,mitob(1),0)
         call brdcst(19+MSGINT,nbf,mitob(1),0)
      endif
      nnp=nbf*(nbf+1)/2
c
c     --- ndmat is the number of density matrices, which is one 
c         less than the number of orbital types
c         virtual orbitals do not contribute. 
      if (mynodeid.eq.0) then
         call iosys('read integer "number of hf density matrices" '//
     $        'from rwf',1,ndmat,0,' ')
      endif
      if(ispar) call brdcst(20+MSGINT,ndmat,mitob(1),0)
c
c     --- what order of derivatives are we to do
      if (mynodeid.eq.0) then
         nderiv=intkey(ops,'nderiv',1,' ')
         if(logkey(ops,'force-constants',.false.,' ')) then
            if(logkey(ops,'force-constants=numerical',.false.,' '))
     $           then
c              numerical force constants, do nothing
            else
c              analytic force constants
               nderiv=2
            endif 
         endif
      endif
      if (ispar) call brdcst(21+MSGINT,nderiv,mitob(1),0)
      if(nderiv.gt.1) then
         call plnkerr('second derivatives not yet implemented',3)
      endif
c
c     --- number of derivatives
      nd1e=3*nat
      nd2e=nd1e*(nd1e+1)/2
c
c     --- allocate some core ---
      f=1
      d=f+ndmat+1
      d2e=d+nnp
      if (nderiv.eq.1) then
         grad=d2e
      else if (nderiv.eq.2) then
         grad=d2e+nd2e
      end if
      dkay=grad+3*nat
      charge=dkay+3*nexch*nbf*nbf
      zan=charge+nat
      c=zan+nat
      top=wpadti(c+3*nat)
c
c     --- retrieve basis set information; returns pointers as well.
      if (mynodeid.eq.0) then
         call iosys('read real "nuclear charges" from rwf',
     $        -1,z(zan),0,' ')
         call iosys('read real coordinates from rwf',-1,z(c),0,' ')
      endif
      if (ispar)then
         call brdcst(22+MSGDBL,z(zan),mdtob(nat),0)
         call brdcst(23+MSGDBL,z(c),mdtob(nat*3),0)
      endif
      call pbasis(nat,nbf,nprim,ncont,ntypes,nbtype,lenxyz,
     $           mxcont,maxl,ptprim,noprim,ptcont,nocont,start,
     $           nocart,nobf,minmom,maxmom,mintyp,nx,ny,nz,
     $           cont,ex,top,z,a)
c
c     --- generate grid points and weights ---
      if (lmax.lt.3 .or. lmax.gt.29)
     $     call plnkerr('m725: invalid lebedev order requested',4)
      nomega=angsiz(lmax)
c
c     --- mxgrd is defaulted to a typical standard grid size.  
      mxgrd=max(mxgrd,(rmax-1)*nomega)
c
c     --- allocate more core.
      ian=top
      ptrad=ian+nat
      radshls=ptrad+rmax*nat
      grid=iadtwp(radshls+nat)
      wts=grid+3*mxgrd
      gradwts=wts+mxgrd
      vwts=gradwts+mxgrd*3*nat
      rnuc=vwts+mxgrd
      amu=rnuc+nat*nat
      pwtx=amu+nat*nat
      rr=pwtx+nat
      radii=rr+nat
      akl=radii+nat
      top=wpadti(akl+nat*nat)
      if (top .gt. canget) then
         if (mynodeid .eq. 0) write(iout,*) 'top,canget',top,canget
         call plnkerr('m511: not enough core to do grid',5)
      endif
c     scratch for dft routines will begin here.
      valuesi=top
      values=iadtwp(valuesi)
      left=iadtwp(canget)-values
c
c     --- finally we know how much we'll have left for kmatrix.
c         for these derivatives, the gradient of the basis functions
c         is always needed.  THIS IS NOT RIGHT
c         increment maxl by 1 for gradients, 2 for laplacian
      bigl=max(maxl+2,2)
c     it needs some arrays that are independent of the grid.
      if(calc.eq.'closed') then
c        (nnp for dtmp, nnp for ktmp, nbf*nbf for kay )
c         dftabs=2*nnp+nbf*nbf
c         those aren't used anymore, but there IS maxl, and there is dsq
         dftabs=nat+nbf*nbf
      else if (calc.eq.'open') then      
c        there are 3*nnp, two for dtmp one for ktmp.  nbf**2 for kay
         dftabs=3*nnp+nbf*nbf
      endif
      kleft=left-dftabs
c
c     --- the remaining core is used to hold arrays which depend
c         on the grid size. we are going to determine the largest
c         grid block we can handle with the available memory to 
c         ensure we can run in the amount of space provided.
c         this means we need to know how much space we need
c         PER GRID POINT. this is given by dftcore
      if(calc.eq.'closed') then
c        dengrida
         dftcore=1
c        fout
         dftcore=dftcore+5 
         if(lyp) then
c           fout needs an extra one
            dftcore=dftcore+1
         endif
c        becke/lyp needs even more
         if(becke.or.lyp) then
c           dengrada,ga
            dftcore=dftcore+3+1
         endif
      else if(calc.eq.'open') then
c        dengrida,dengridb
         dftcore=2
c        fout
         dftcore=dftcore+5 
         if(lyp) then
            dftcore=dftcore+1
         endif
c        becke/lyp needs even more
         if(becke.or.lyp) then
c           dengrada,dengradb,ga,gb,queue
            dftcore=dftcore+3+3+1+1+3
            if(lyp) then
c              gab
               dftcore=dftcore+1
            endif
         endif
      endif
c
c     --- both open/closed need
c         phi,grad,hess,scr,tea,tmpgwt,nzptrs
c           (treat this as if were a real*8 array for now)
c         accum
      dftcore=dftcore+nbf+3*nbf+6*nbf+1+nbf+1+1+nbf*3
c
c     --- this takes us up to the beginning of itch in kmatrix.
c         the functionals (except slater) need some scratch space
      dfttmp=0
      if (becke) dfttmp=5
      if (vwn) then
         dfttmp=max(dfttmp,16)
      else if (lyp) then
         dfttmp=max(dfttmp,12)
      endif
c
c     --- direct k need some core too
c         needs rsq,s,r,tee,xyzpow
      dfttmp=max(dfttmp,1+3*mxcont+3*bigl)
c
c     ta da..... this is how much room we need PER GRID POINT
      dftcore=dftcore+dfttmp
c
c     --- finally, determine the maximum block size which will fit
      mxgbsiz=kleft/dftcore
c
c     let user override, to avoid using all the memory there is and 
c     choking the machine, but only to DECREASE it, not to increase it.
      if (mynodeid .eq. 0) then
         mxgbsiz=min(intkey(ops,'scf=mxgbsiz',mxgbsiz,' '),mxgbsiz)
         mxgbsiz=min(mxgrd,mxgbsiz)
      endif
      if (ispar) call brdcst(24+MSGINT,mxgbsiz,mitob(1),0)
      if(prnt .and. (mynodeid .eq.0)) then
         write(iout,1070) mxgbsiz
      endif
c
c     --- pick up the atomic numbers and density matrices
      if (mynodeid.eq.0) then
         call iosys('read integer "atomic numbers" from rwf',
     $        -1,a(ian),0,' ')
         call iosys('read real f from rwf',ndmat+1,z(f),0,' ')
         call iosys('read real "hf density matrix" from rwf ',
     $        nnp,z(d),0,' ')
      endif
      if (ispar) then
         call brdcst(25+MSGINT,a(ian),mitob(nat),0)
         call brdcst(26+MSGDBL,z(f),mdtob(ndmat+1),0)
         call brdcst(27+MSGDBL,z(d),mdtob(nnp),0)
      endif
c
c     --- calculate derivative exchange-correlation integrals ---
c         and the contribution to the Exc gradient
      call dxcint(z(values),z(d),nbf,nnp,z(dkay),nexch,z(wts),
     $            ndmat,nat,mxgrd,exc,slater,becke,lyp,vwn,
     $            calce,calc,dmcut,dencut,eexch,prteexch,
     $            mxgblk,mxgbsiz,a(valuesi),
     $            z(c),z(ex),z(cont),a(ptprim),a(noprim),
     $            a(nocont),a(ptcont),mxcont,nprim,
     $            ntypes,nbtype,ncont,a(start),
     $            a(nocart),a(nobf),a(maxmom),a(minmom),a(mintyp),
     $            a(nx),a(ny),a(nz),z(grid),z(charge),bigl,
     $            ops,nderiv,z(grad),z(gradwts),a(ian),a(ptrad),
     $            a(radshls),z(vwts),z(rnuc),z(amu),z(pwtx),z(rr),
     $            z(radii),z(akl),rmax,lmax,nomega,nradial,grdtyp,
     $            adjust,minesz)
c
c     --- check the energies ---
      if (mynodeid .eq. 0) then
         call iosys('read real "hf xc energy" from rwf',-1,exc511,0,' ')
         if(abs((exc511-exc)/exc511).gt.1.0d-06) then
            write (iout,87) exc,exc511
 87         format (/5x,'calculated xc energy is:',g20.12,
     $           /,5x,'         xc from scf is:',g20.12)
c            call plnkerr('energies do not agree !!!',6)
         end if
c
c     --- store these contributions on the rwf ---
         call iosys('write real "xc integral derivatives" to rwf',
     $        3*nat,z(grad),0,' ')
         if(nderiv.eq.2) then
            call iosys('write real "xc integral force constants"'
     $           //' to rwf',nd2e,z(d2e),0,' ')
         endif
c
c     --- print gradient contributions and form total gradient ---
         if (logkey(ops,'print=gradient=xc',.false.,' ')) then
            write (iout,1010)
            call matout(z(grad),3,nat,3,nat,iout)
            if (nderiv.eq.2) then
               write (iout,1020)
               call print(z(d2e),nd2e,nd1e,iout)
            end if
         end if
c
c     --- form total two-electron gradient ---
         call iosys(
     $        'read real "coulomb integral derivatives" from rwf',
     $        3*nat,z(zan),0,' ')
         call vadd(z(grad),z(grad),z(zan),3*nat)
         if(nderiv.eq.2) then
            call iosys('read "coulomb integral force constants"'
     $           //' from rwf',nd2e,z(zan),0,' ')
            call vadd(z(d2e),z(d2e),z(zan),nd2e)
         endif
c     
         if (logkey(ops,'print=gradient=two-electron',.false.,' '))
     $        then
            write (iout,1030)
            call matout(z(grad),3,nat,3,nat,iout)
            if (nderiv.ge.2) then
               write (iout,1040)
               call print(z(d2e),nd2e,nd1e,iout)
            end if
         end if
c
c     --- form total gradients ---
         call iosys('read real "one-electron derivatives" from rwf',
     $        3*nat,z(zan),0,' ')
         call vadd(z(grad),z(grad),z(zan),3*nat)
         if (nderiv.ge.2) then
            call iosys(
     $           'read real "one-electron force-constants" from rwf',
     $           nd2e,z(zan),0,' ')
            call vadd(z(d2e),z(d2e),z(zan),nd2e)
         end if
c
         if (logkey(ops,'print=gradient=total',.false.,' ')) then
            write (iout,1050)
            call matout(z(grad),3,nat,3,nat,iout)
            if (nderiv.ge.2) then
               write (iout,1060)
               call print(z(d2e),nd2e,nd1e,iout)
            end if
         end if
c
c     --- save the gradients to the rwf ---
         call iosys('write real "hf first derivatives" to rwf',
     $        3*nat,z(grad),0,' ')
         call iosys('write real "cartesian first derivatives" to rwf',
     $        3*nat,z(grad),0,' ')
         if (nderiv.ge.2) then
            call iosys('write real "integral force constants" to rwf',
     $           nd2e,z(d2e),0,' ')
         end if
c
c     --- and exit gracefully ---
         call chainx(0)
      endif
c
c
      return
      end
