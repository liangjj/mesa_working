*deck @(#)vibfrc.f	5.1 11/6/94
      subroutine vibfrc(fsq,natoms,ian,atmass,cmass,
     $                  nvar,cref,natoms3,nnp,zf,vname)
c***begin prologue     vibfrc.f
c***date written       yymmdd  
c***revision date      11/6/94      
c
c***keywords           
c***author             page, michael(nrl)
c***source             @(#)vibfrc.f	5.1   11/6/94
c***purpose            
c***description
c
c          mass weights the cartesian force constant matrix, projects
c          rotations and translations out and diagonalizes it, yielding
c          vibrational frequencies.
c***references
c
c***routines called
c
c***end prologue       vibfrc.f
      implicit none
c     --- input variables -----
      integer natoms,nvar,natoms3,nnp
c     --- input arrays (unmodified) ---
      integer ian(natoms)
      character*(*) vname(nvar)
      real*8 atmass(natoms)
c     --- input arrays (scratch) ---
c     --- output arrays ---
c     --- output variables ---
c     --- scratch arrays ---
      real*8 cmass(natoms3)
      real*8 fsq(natoms3,natoms3)
      real*8 cref(natoms3)
      real*8 zf(*)
c     --- local variables ---
      integer d2ecycl
      integer inp,iout
      integer fctri,eigval,eigvec,t1,t2,t3,projop,proj
      integer t4,t5,t6,top,maxcor,jatom
      integer k,icfx,b
      integer wpadti
      logical prnt,chkpt,singpt,cartfx
      logical projf
      logical debug
      real*8 energy,rmax,rmin,rlim,stpsize
      real*8 zpe
c
      parameter (debug=.false.)
c
      common/d2einf/energy,rmax,rmin,rlim,d2ecycl,
     $               prnt,chkpt,singpt,stpsize,cartfx
      common /io/ inp,iout
c
 1010 format('  vibfrc: top = ',i8,'maxcor=',i8)
 1020 format(5x,'cartesian force constant matrix')
 1040 format(5x,'mass weighted fx before freq1 ')
 1050 format (5x,'unprojected frequencies')
 1060 format(//,' mass weighted fx after freq1 ',/)
 1070 format(/' vibfrc before proj : masses ',/,3f18.10)
 1080 format(5x,
     $   'projected mass weighted cartesian force constant matrix ')
 1090 format (5x,'frequencies with rotations and translations ',
     $           'projected out')
 1100 format(5x,'zero point vibrational energy: ',f8.2,' kcal/mole')
 1110 format(5x,'***** normal modes of vibration ***** ')
c
c     --- allocate some core
      fctri=1
      eigval=fctri+nnp
      eigvec=eigval+natoms3
      t1=eigvec+natoms3**2
      t2=t1+natoms3**2
      t3=t2+natoms3**2
      projop=t3+nnp
      proj=projop+6*natoms3
      t4=proj+natoms3**2
      t5=t4+natoms3**2
      t6=t5+natoms3**2
      b=t6+nvar**2
      top=wpadti(b+nvar*natoms3)
c
      call getscm(top,zf,maxcor,'m205',0)
      if(top.ge.maxcor) then
         write(iout,1010) top,maxcor
         call lnkerr(' need more space ')
      endif
c
c     --- fill up the cartesian mass array
      do 10 jatom=1,natoms
         cmass((jatom-1)*3+1)=atmass(jatom)
         cmass((jatom-1)*3+2)=atmass(jatom)
         cmass((jatom-1)*3+3)=atmass(jatom)
   10 continue
c
c    --- put force constant matrix into a triangle and then
c        write it to rwf.  then overwrite the triangle with
c        the mass weighted force constant matrix.
      call sqtotr(zf(fctri),fsq,natoms3,nnp)
      call iosys('write real "cartesian second derivatives" to rwf',
     $            nnp,zf(fctri),0,' ')
      if(debug) then
         write(iout,1020)
         call matout(fsq,natoms3,natoms3,natoms3,natoms3,iout)
      endif
c
      call masswt(fsq,cmass,natoms3)
      call sqtotr(zf(fctri),fsq,natoms3,nnp)
c
      if(debug) then
         write(iout,1040)
         call print(zf(fctri),nnp,natoms3,iout)
         write (iout,1050)
      endif
      call freq(zf(fctri),zf(eigval),zf(eigvec),zf(t1),zf(t2),
     $          zf(t3),natoms3,nnp,zpe)
      if(debug) then
         write(iout,1110)
         call matout(zf(eigvec),natoms3,natoms3,natoms3,natoms3,iout)
         write(iout,1060)
         call print(zf(fctri),nnp,natoms3,iout)
      endif
c
c
      projf=.true.
      if(projf) then
         if(debug) then
            write(iout,1070) (atmass(k),k=1,natoms)
         endif
         call project(zf(fctri),zf(t3),zf(projop),zf(proj),cref,
     $                atmass,zf(t4),zf(t5),natoms,natoms3,nnp)
         if(debug) then
            write(iout,1080)
            call print(zf(t3),nnp,natoms3,iout)
            write (iout,1090)
         endif
         call freq(zf(t3),zf(eigval),zf(eigvec),zf(t1),zf(t2),
     $             zf(fctri),natoms3,nnp,zpe)
      endif
c
c     --- print the zero point energy
      write(iout,1100) zpe
      write(iout,1110)
c
      icfx=1
      call eigint(zf(eigvec),zf(t2),natoms3,nvar,vname,zf(b),
     $            zf(t6),icfx)
c
c
      return
      end
