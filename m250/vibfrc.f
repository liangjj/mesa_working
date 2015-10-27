*deck @(#)vibfrc.f	5.1  11/6/94
      subroutine vibfrc(fsq,ftri,natoms,ian,atmass,cmass,
     $             cref,natoms3,nnp,zf,order,vname,nvar,
     $             projg,sstep,phase,nbf,saddle)
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
c   mass weights the cartesian force constant matrix, projects out
c   rotations translations, and possibly the gradient, diagonalizes
c   it, yielding vibrational frequencies.
c
c***references
c
c***routines called
c
c***end prologue       vibfrc.f
      implicit none
c     --- input variables -----
      integer natoms,natoms3,nnp,nvar,nbf
      real*8 sstep,phase
c     --- input arrays (unmodified) ---
      integer ian(natoms),order(natoms3)
      character*(*) vname(nvar)
c     --- input arrays (scratch) ---
      real*8 atmass(natoms),cmass(natoms3)
      real*8 fsq(natoms3,natoms3),ftri(nnp)
      real*8 cref(natoms3)
      real*8 zf(*)
c     --- output arrays ---
c     --- output variables ---
c     --- scratch arrays ---
c     --- local variables ---
      integer inp,iout
      integer eigval,eigvec,trvec,t1,t2,t3,projop,proj
      integer t4,t5,binv,scr6,scr7,b,xstep,cmo,cphf,top
      integer jatom,icfx,gcm
      real*8 zpe
      logical projg,saddle,debug
c
      parameter (debug=.true.)
c
      common /io/ inp,iout
 1010 format('  vibfrc: top = ',i8)
 1020 format(/' cartesian force constant matrix'//)
 1030 format(/' vibfrc: masses ',3f18.10)
 1040 format(//,' mass weighted fx before freq1 ',/)
 1050 format (//,'      frequencies with rotations and translations '
     $   ,       'projected out',//)
 1060 format (//,'      frequencies with rotations, translations '
     $           ,       'and gradient projected out',//)
 1070 format(/' zero point vibrational energy: ',f8.2,' kcal/mole')
 1080 format(/'  *** normal modes of vibration *** '/)
c
c     -----allocate some core
      eigval=1
      eigvec=eigval+natoms3
      trvec=eigvec+natoms3**2
      t1=trvec+natoms3
      t2=t1+natoms3**2
      t3=t2+natoms3**2
      projop=t3+nnp
      proj=projop+7*natoms3
      t4=proj+natoms3**2
      t5=t4+natoms3**2
      binv=t5+natoms3**2
      scr6=binv+natoms3*nvar
      scr7=scr6+2*nvar
      b=scr7+nvar*nvar
      xstep=b+nvar*3*natoms
      cmo=xstep+natoms3
      cphf=cmo+nbf**2
      top=cphf+nbf**2
      if(top.ge.10000) then
         write(iout,1010) top
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
      if(saddle) then
         call trtosq(fsq,ftri,natoms3,nnp)
         write(iout,1020)
         call matout(fsq,natoms3,natoms3,natoms3,natoms3,iout)
c
         call masswt(fsq,cmass,natoms3)
         call sqtotr(ftri,fsq,natoms3,nnp)
c
         if(debug) then
            write(iout,1040)
            call print(ftri,nnp,natoms3,iout)
         endif
         call project(ftri,zf(t3),zf(projop),zf(proj),cref,
     $                atmass,zf(t4),zf(t5),natoms,natoms3,nnp,
     $                projg)
c
         if(projg) then
            write (iout,1060)
         else
            write (iout,1050)
         endif
c
         call freq(zf(t3),zf(eigval),zf(eigvec),zf(t1),zf(t2),
     $             ftri,natoms3,nnp,zpe)
         write(iout,1070) zpe
         write(iout,1080)
c        --- dont calculate internal coordinate normal modes
c            because the ctoz transformation matrix will be
c            wrong during a path following
         icfx=1
         call eigint(zf(eigvec),zf(t2),natoms3,nvar,vname,zf(b),
     $               zf(scr7),icfx,zf(binv),order,zf(scr6),
     $               zf(trvec))
      end if
c
      gcm=projop+6*natoms3
      if(.not.saddle) then
         call getgrd(zf(gcm),atmass,natoms,natoms3)
      endif
c
      call stepp(zf(gcm),zf(xstep),cref,zf(trvec),atmass,
     $           natoms,natoms3,projg,sstep,phase)
c
c
      return
      end
