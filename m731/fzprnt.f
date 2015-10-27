*deck @(#)fzprnt.f	5.1  11/6/94
      subroutine fzprnt(nz,ianz,iz,f)
c***begin prologue     fzprnt.f
c***date written       850601  yymmdd
c***revision date      11/6/94
c***keywords           z-matrix, forces
c***author             binkley, et al., gaussian82
c                      martin, richard (lanl)
c***source             @(#)fzprnt.f	5.1   11/6/94
c***purpose            prints the forces in a z-matrix format.
c***description
c
c     call fzprnt(nz,ianz,iz,f)
c
c
c     routine to print the internal coordinate forces in a
c     z-matrix like manner.  this routine produces output very
c     similar to zprint, but differs in the following:
c       1.  formats produce more significant figures.
c       2.  fzprnt knows how to get forces from a single linear array.
c
c     arguments:
c
c       nz     ... number of lines in the z-matrix.
c       ianz   ... integer array of length nz containing the atomic numbers
c                  of the z-matrix centers.
c       iz     ... the integer connectivity information associated
c                  with the z-matrix.
c       f      ... real array containing the internal-coordinate forces,
c                  stored in the same arrangement that variables are
c                  numbered in a z-matrix.
c
c***references
c***routines called    fillel(util)
c***end prologue       fzprnt.f
      implicit none
c     --- input variables -----
      integer nz
c     --- input arrays (unmodified) ---
      integer ianz(nz),iz(4,nz)
      real*8 f(nz)
c     --- input arrays (scratch) ---
c     --- output arrays ---
c     --- output variables ---
c     --- scratch arrays ---
c     --- local variables ---
      integer inp,iout
      integer icard,np1,np2,np3,idx,icent
      integer maxel
      character el(106)*2
      data maxel/104/, el(1)/'x'/
      save maxel,el
c
      common/io/inp,iout
c
 1010 format(5x,'internal coordinate forces (hartrees/bohr,radian)')
 1020 format(5x,'   cent atom   n1',5x,'length',8x,'n2',5x,'alpha',
     $       8x,'n3',7x,'beta',7x,'j')
 2110 format(5x,3x,2x,i2,3x,a2)
 2120 format(5x,3x,      7x,a2)
 2210 format(5x,3x,2x,i2,3x,a2,2x,i3,1x,f10.6,' (',i3,')')
 2220 format(5x,3x,      7x,a2,2x,i3,1x,f10.6,' (',i3,')')
 2310 format(5x,3x,2x,i2,3x,a2,2x,i3,1x,f10.6,' (',i3,') ',
     $          1x,i2,1x,f10.6,' (',i3,')')
 2320 format(5x,3x,      7x,a2,2x,i3,1x,f10.6,' (',i3,') ',
     $          1x,i2,1x,f10.6,' (',i3,')')
 2410 format(5x,3x,2x,i2,3x,a2,2x,i3,1x,f10.6,' (',i3,') ',
     $          1x,i2,1x,f10.6,' (',i3,') ',i2,1x,f10.6,' (',i3,') ',
     $          i2)
 2420 format(5x,3x,      7x,a2,2x,i3,2x,f10.6,' (',i3,') ',
     $          1x,i2,1x,f10.6,' (',i3,') ',i2,1x,f10.6,' (',i3,') ',
     $          i2)
c
c     --- print the heading.
      write(iout,1010)
      write(iout,1020)
c
c     --- first card of the z-matrix.
      if(nz.ge.1) then
         call fillel(0,maxel,el(2))
         idx=ianz(1)+2
         if(ianz(1).ge.0) then
            icent=1
            write(iout,2110) icent,el(idx)
         else
            icent=0
            write(iout,2120) el(idx)
         endif
      endif
c
c     --- second card of the z-matrix.
      if(nz.ge.2) then
         np1=1
         idx=ianz(2)+2
         if(ianz(2).ge.0) then
            icent=icent+1
            write(iout,2210) icent,el(idx),iz(1,2),f(np1),np1
         else
            write(iout,2220)       el(idx),iz(1,2),f(np1),np1
         endif
      endif
c
c     --- third card.
      if(nz.ge.3) then
         np1=2
         np2=nz
         idx=ianz(3)+2
         if(ianz(3).ge.0) then
            icent=icent+1
            write(iout,2310) icent,el(idx),iz(1,3),f(np1),np1,
     $                       iz(2,3),f(np2),np2
         else
            write(iout,2320)       el(idx),iz(1,3),f(np1),np1,
     $                       iz(2,3),f(np2),np2
         endif
      endif
c
c     --- cards 4 through nz.
      if(nz.ge.4) then
         do 100 icard=4,nz
            np1=icard-1
            np2=nz+icard-3
            np3=nz*2+icard-6
            idx=ianz(icard)+2
            if(ianz(icard).ge.0) then
               icent=icent+1
               write(iout,2410) icent,el(idx),iz(1,icard),
     $                          f(np1),np1,iz(2,icard),f(np2),
     $                          np2,iz(3,icard),f(np3),np3,
     $                          iz(4,icard)
            else
               write(iout,2420)       el(idx),iz(1,icard),
     $                          f(np1),np1,iz(2,icard),f(np2),
     $                          np2,iz(3,icard),f(np3),np3,
     $                          iz(4,icard)
            endif
  100    continue
      endif
c
c
      return
      end
