*deck @(#)frmain.f	5.1  11/6/94
      subroutine frmain(nvar,nz,natoms,ops,vname,ian,atmass,
     $                  cmass,xc,zf,ftric,fsq,nvvc,order)
c***begin prologue     frmain.f
c***date written       yymmdd  
c***revision date      11/6/94      
c
c***keywords           
c***author             page, michael(nrl)
c***source             @(#)frmain.f	5.1   11/6/94
c***purpose            
c***description
c
c     local arrays of note.
c     vname  ... variable names vector; character*16(nvar).
c     ian    ... atomic numbers
c     atmass ... atomic masses
c
c***references
c
c***routines called
c
c***end prologue       frmain.f
      implicit none
c     --- input variables -----
      integer nvar,nz,natoms,nvvc
c     --- input arrays (unmodified) ---
      character*(*) ops,vname(nvar)
      integer ian(nz),order(3*natoms)
c     --- input arrays (scratch) ---
      real*8 atmass(natoms)
      real*8 ftric(nvvc),cmass(3*natoms)
      real*8 fsq(3*natoms,3*natoms)
      real*8 zf(2),xc(3*natoms)
c     --- output arrays ---
c     --- output variables ---
c     --- scratch arrays ---
c     --- local variables ---
      integer inp,iout
      integer nat3,nwords
c
      common/io/inp,iout
c
 1000 format(1x,'m204:calculation of vibrational frequencies')
c
c     --- announce our presence.
      write(iout,1000)
c
      nat3=natoms*3
      call iosys('read real coordinates from rwf',nat3,xc,0,' ')
c     --- get force constants from chk file
c     call readfx(ftric,nvvc)
      call iosys('read real "cartesian second derivatives" from rwf',
     $            nvvc,ftric,0,' ')
      call iosys('read integer "atomic numbers" from rwf',
     $           -1,ian,0,' ')
      call filmas(0,iout,ian,natoms,.false.,atmass,' ')
c
c     --- read in the variable names
      nwords=nvar*len(vname(1))
      call iosys('read character "variable names" from rwf',
     $            nwords,0,0,vname)
c
c     --- compute the frequencies
      call vibfrc(fsq,ftric,natoms,ian,atmass,cmass,xc,
     $                   nat3,nvvc,zf,order,vname,nvar)
c
c
      return
      end
