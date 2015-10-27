*deck  @(#)frmain.f	5.1 11/6/94
      subroutine frmain(nvar,nz,natoms,ops,vname,ian,atmass,
     $                  cmass,xc,zf,ftric,fsq,nvvc,order,projg,
     $                  sstep,phase,saddle)
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
      logical saddle,projg
      real*8 sstep,phase
c     --- input arrays (unmodified) ---
      integer ian(nz),order(3*natoms)
      character*(*) ops,vname(nvar)
c     --- input arrays (scratch) ---
      real*8 atmass(natoms),ftric(nvvc),cmass(3*natoms)
      real*8 fsq(3*natoms,3*natoms)
      real*8 zf(*),xc(3*natoms)
c     --- output arrays ---
c     --- output variables ---
c     --- scratch arrays ---
c     --- local variables ---
      integer inp,iout
      integer nat3,nwords,nbf
c
      common/io/inp,iout
c
 2000 format(/,5x,'calculation of vibrational frequencies'/)
c
c     --- get current coordinates
      nat3=natoms*3
      call iosys('read real coordinates from rwf',nat3,xc,0,' ')
c
c     --- perhaps pick up second derivatives.
      if(saddle) then
         call iosys('read real "cartesian second derivatives" from rwf',
     $               nvvc,ftric,0,' ')
      end if
      call iosys('read integer "atomic numbers" from rwf',
     $          -1,ian,0,' ')
      call filmas(0,iout,ian,natoms,.false.,atmass,' ')
c
c     --- read in the variable names,etc.
      nwords=nvar*len(vname(1))
      call iosys('read character "variable names" from rwf',
     $            nwords,0,0,vname)
      call iosys('read integer "number of basis functions" from rwf',
     $           1,nbf,0,' ')
c
c     --- get frequencies
      call vibfrc(fsq,ftric,natoms,ian,atmass,cmass,xc,
     $            nat3,nvvc,zf,order,vname,nvar,projg,
     $            sstep,phase,nbf,saddle)
c
c
      return
      end
