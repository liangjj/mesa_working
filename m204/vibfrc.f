*deck @(#)vibfrc.f	5.1 11/6/94
      subroutine vibfrc(fsq,ftri,natoms,ian,atmass,cmass,
     $                  cref,natoms3,nnp,zf,order,vname,nvar)
c***begin prologue     vibfrc.f
c***date written       yymmdd  
c***revision date      11/6/94      
c
c***keywords           
c***author             page, michael(nrl)
c***source             @(#)vibfrc.f	5.1   11/6/94
c***purpose            
c***description
c     mass weights the cartesian force constant matrix, projects
c     rotations and translations out and diagonalizes it, yielding
c     vibrational frequencies.
c
c***references
c
c***routines called
c
c***end prologue       vibfrc.f
      implicit none
c     --- input variables -----
      integer natoms,natoms3,nnp,nvar
c     --- input arrays (unmodified) ---
      integer ian(natoms),order(natoms3)
      character*(*) vname(nvar)
      real*8 ftri(nnp),atmass(natoms)
c     --- input arrays (scratch) ---
c     --- output arrays ---
      real*8 cmass(natoms3)
c     --- output variables ---
c     --- scratch arrays ---
      real*8 fsq(natoms3,natoms3)
      real*8 cref(natoms3)
      real*8 zf(*)
c     --- local variables ---
      integer inp,iout
      integer eigval,eigvec,t1,t2,t3,projop,proj,t4,t5
      integer binv,scr6,scr7,b,top
      integer wpadti
      integer jatom,ndegf,icfx,maxcor,k
      logical projf
      logical debug
      real*8 zpe
c
      parameter (debug=.false.)
c
      common /io/ inp,iout
c
 1010 format('  vibfrc: top = ',i8,'maxcor=',i8)
 1020 format(5x,'cartesian force constant matrix')
 1030 format(5x,'vibfrc:masses ',3f18.10)
 1040 format(5x,'mass weighted fx before freq1 ')
 1050 format (5x,'unprojected frequencies')
 1060 format(5x,'***** normal modes of vibration ***** ')
 1070 format(5x,'mass weighted fx after freq1 ')
 1080 format(5x,'vibfrc before proj:masses ',/,3f18.10)
 1090 format(5x,
     $       'projected mass weighted cartesian force constant matrix')
 1100 format (5x,'frequencies with rotations and translations ',
     $          'projected out')
 1110 format(5x,'zero point vibrational energy: ',f8.2,' kcal/mole')
c
c     --- allocate some core
      eigval=1
      eigvec=eigval+natoms3
      t1=eigvec+natoms3**2
      t2=t1+natoms3**2
      t3=t2+natoms3**2
      projop=t3+nnp
      proj=projop+6*natoms3
      t4=proj+natoms3**2
      t5=t4+natoms3**2
      binv=t5+natoms3**2
      scr6=binv+natoms3*nvar
      scr7=scr6+2*nvar
      b=scr7+nvar*nvar
      top=wpadti(b+nvar*3*natoms)
c
      call getscm(top,zf,maxcor,'m204',0)
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
      call trtosq(fsq,ftri,natoms3,nnp)
      write(iout,1020)
      call matout(fsq,natoms3,natoms3,natoms3,natoms3,iout)
c
      if(debug) then
         write(iout,1030) (atmass(k),k=1,natoms)
      endif
      call masswt(fsq,cmass,natoms3)
      call sqtotr(ftri,fsq,natoms3,nnp)
      if(debug) then
         write(iout,1040)
         call print(ftri,nnp,natoms3,iout)
      endif
c
      if(debug) then
         write (iout,1050)
      endif
      call freq(ftri,zf(eigval),zf(eigvec),zf(t1),zf(t2),
     $          zf(t3),natoms3,nnp,zpe)
      if(debug) then
         write(iout,1060)
         call matout(zf(eigvec),natoms3,natoms3,natoms3,natoms3,iout)
         write(iout,1070)
         call print(ftri,nnp,natoms3,iout)
      endif
c
c
      projf=.true.
      if(projf) then
         if(debug) then
            write(iout,1080) (atmass(k),k=1,natoms)
         endif
         call project(ftri,zf(t3),zf(projop),zf(proj),cref,
     $                atmass,zf(t4),zf(t5),natoms,natoms3,nnp)
c
         if(debug) then
            write(iout,1090)
            call print(zf(t3),nnp,natoms3,iout)
            write (iout,1100)
         endif
c
         call freq(zf(t3),zf(eigval),zf(eigvec),zf(t1),zf(t2),
     $             ftri,natoms3,nnp,zpe)
      endif
c
c
      write(iout,1110) zpe
      write(iout,1060)
c
c     --- if a full set of internal coordinates is available, then
c         transform normal modes to internals too.
      ndegf=natoms3-6
      if(nvar.ge.ndegf) then
         icfx=0
      else
         icfx=1
      endif
      call eigint(zf(eigvec),zf(t2),natoms3,nvar,vname,zf(b),
     $            zf(scr7),icfx,zf(binv),order,zf(scr6))
c
c
      return
      end
