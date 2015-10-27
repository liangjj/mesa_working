*deck @(#)vibfreq.f	5.1  11/6/94
      subroutine vibfreq(frcnst,natoms,ian,atmass,cmass,
     $                   xx,ff,nvv,nvar,maxpt,b,cref,natoms3,nnp,
     $                   iscr11,zf,vname)
c***begin prologue     vibfreq.f
c***date written       yymmdd  
c***revision date      11/6/94      
c
c***keywords           
c***author             
c***source             @(#)vibfreq.f	5.1   11/6/94
c***purpose            
c***description
c      this routine:
c      1.  transforms the internal coordinate force constant matrix to
c          cartesian coordinates. this generates a 3n-6 x 3n-6 cartesian
c          force constant matrix. the six rows and columns remaining are
c          those cartesian directions held fixed in m202 (see subroutine
c          ztoc). these are found by considering the rotational and
c          translational invariance of the energy.
c      2.  mass weights the cartesian force constant matrix, projects
c          rotations and translations out and diagonalizes it, yielding
c          vibrational frequencies.
c
c***references
c
c***routines called
c
c***end prologue       vibfreq.f
      implicit none
      integer natoms,nvv,nvar,maxpt,natoms3,nnp
c     --- input variables -----
c     --- input arrays (unmodified) ---
      integer ian(natoms)
      character *(*) vname(nvar)
      real*8 atmass(natoms),frcnst(nvv)
      real*8 b(natoms3,nvar),cref(natoms3)
c     --- input arrays (scratch) ---
c     --- output arrays ---
c     --- output variables ---
c     --- scratch arrays ---
      integer iscr11(natoms3)
      real*8 xx(nvar,maxpt),cmass(natoms3)
      real*8 ff(nvar,maxpt)
      real*8 zf(*)
c     --- local variables ---
      integer inp,iout
      integer d2ecycl
      integer fsq,fc,tlxb,fctri,eigval,eigvec,t1,t2,t3
      integer projop,proj,t4,t5,bvec,l,t6,t7,t8,t9
      integer t10,t11,t12,t13,top
      integer maxcor,wpadti
      integer i,j,k,kk,jatom,icfx
      logical prnt,chkpt,singpt,cartfx
      logical debug
      logical projf
      real*8 energy,rmax,rmin,rlim,stpsize
      real*8 zpe
c
      parameter (debug=.false.)
c
      common/d2einf/energy,rmax,rmin,rlim,d2ecycl,
     $               prnt,chkpt,singpt,stpsize,cartfx
      common /io/ inp,iout
c
 1000 format('  vibfreq: top = ',i8)
 1010 format(5x,'reference cartesian geometry')
 1020 format(10x,3f18.5)
 1030 format(5x,'b matrix')
 1040 format(/' unfilled force constant matrix'//)
 1050 format(' coordinates before fillfc ',(3f20.8))
 1060 format(' calling fillfc zf addresses',/,(i20))
 1070 format(5x,'cartesian force constant matrix')
 1080 format(5x,'mass weighted fx before freq1')
 1090 format (5x,'unprojected frequencies')
 1100 format(//,' mass weighted fx after freq1 ',/)
 1110 format(5x,
     $  'projected mass weighted cartesian force constant matrix ')
 1120 format (5x,'frequencies with rotations and translations ',
     $           'projected out')
 1130 format(5x,'zero point vibrational energy: ',f8.2,' kcal/mole')
 1140 format(5x,'***** normal modes of vibration *****')
c
c     --- allocate some core
      fsq=1
      fc=fsq+nvar*nvar
      tlxb=fc+natoms3**2
      fctri=tlxb+natoms3*nvar
      eigval=fctri+nnp
      eigvec=eigval+natoms3
      t1=eigvec+natoms3**2
      t2=t1+natoms3**2
      t3=t2+natoms3**2
      projop=t3+nnp
      proj=projop+6*natoms3
      t4=proj+natoms3**2
      t5=t4+natoms3**2
      bvec=t5+natoms3**2
      l=bvec+6*natoms3
      t6=l+natoms3*6
      t7=t6+6
      t8=t7+6*natoms3
      t9=t8+nvar**2
      t10=t9+natoms3
      t11=t10+6
      t12=t11+6*natoms3
      t13=t12+6*natoms3
      top=wpadti(t13+6*natoms3)
c
      call getscm(top,zf,maxcor,'m205',0)
      if(top.ge.maxcor) then
         write(iout,1000) top
         call lnkerr(' need more space ')
      endif
c
c     --- reference geometry
      if(debug) then
         write(iout,1010)
         do 44 i=1,natoms
            k=3*(i-1)+1
            kk=k+2
            write(iout,1020) (cref(j),j=k,kk)
   44    continue
      endif
c
c     --- fill up the cartesian mass array
      do 10 jatom=1,natoms
         cmass((jatom-1)*3+1)=atmass(jatom)
         cmass((jatom-1)*3+2)=atmass(jatom)
         cmass((jatom-1)*3+3)=atmass(jatom)
   10 continue
c
c     --- dump the b matrix
      if(debug) then
         write(iout,1030)
         call matout(b,natoms3,nvar,natoms3,nvar,iout)
      endif
c
c
      call trtosq(zf(fsq),frcnst,nvar,nvv)
      call ebct(zf(tlxb),zf(fsq),b,nvar,nvar,natoms3)
      call ebc(zf(fc),b,zf(tlxb),natoms3,nvar,natoms3)
c
c     --- fill out the cartesian force constant matrix, but first
c         make sure that the appropriate cartesian directions were
c         fixed during the finite difference procedure. that is:
c         x,y and z on atom1; x and y on atom 2 and z on atom3
c         these rows and columns of the transformed force constant
c         matrix should be zero
      if(debug) then
         write(iout,1040)
         call matout(zf(fc),natoms3,natoms3,natoms3,natoms3,iout)
         write(iout,1050) (cref(i),i=1,natoms3)
      endif
c
c     --- transform the gradient
      call embc(zf(t9),b,ff,natoms3,nvar,1)
      if(debug) then
         write(iout,1060) fc,l,bvec,t4,t5,t6,t7,t9,t10,t11,t12,t13
      endif
      call fillfc(zf(fc),zf(l),natoms3,natoms,
     $            zf(bvec),ff,cref,nvar,maxpt,
     $            zf(t4),zf(t5),iscr11,zf(t6),zf(t7),
     $            zf(t9),zf(t10),zf(t11),zf(t12),zf(t13))
      if(debug) then
         write(iout,1070)
         call matout(zf(fc),natoms3,natoms3,natoms3,natoms3,iout)
      endif
c
c     --- put force constant matrix into a triangle and then
c         write it to rwf.  then overwrite the triangle with
c         the mass weighted force constant matrix.
      call sqtotr(zf(fctri),zf(fc),natoms3,nnp)
      call iosys('write real "cartesian second derivatives" to rwf',
     $            nnp,zf(fctri),0,' ')
      call masswt(zf(fc),cmass,natoms3)
      call sqtotr(zf(fctri),zf(fc),natoms3,nnp)
      if(debug) then
         write(iout,1080)
         call print(zf(fctri),nnp,natoms3,iout)
         write (iout,1090)
         call freq(zf(fctri),zf(eigval),zf(eigvec),zf(t1),zf(t2),
     $             zf(t3),natoms3,nnp,zpe)
         write(iout,1100)
         call print(zf(fctri),nnp,natoms3,iout)
      endif
c
c     --- project out translations/rotations
      projf=.true.
      if(projf) then
         call project(zf(fctri),zf(t3),zf(projop),zf(proj),cref,
     $                atmass,zf(t4),zf(t5),natoms,natoms3,nnp)
         if(debug) then
            write(iout,1110)
            call print(zf(t3),nnp,natoms3,iout)
c           write (iout,1120)
         endif
         call freq(zf(t3),zf(eigval),zf(eigvec),zf(t1),zf(t2),
     $             zf(fctri),natoms3,nnp,zpe)
      endif
c
c     --- print the zero point energy
      write(iout,1130) zpe
      write(iout,1140)
c
c     --- print normal modes
      icfx=0
      call eigint(zf(eigvec),zf(t2),natoms3,nvar,vname,b,
     $            zf(t8),icfx)
c
c
      return
      end
