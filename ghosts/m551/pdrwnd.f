*deck %W%  %G%
      subroutine pdrwnd ( nfti )
C
C***Begin prologue
C***Date written       871022   (yymmdd)
C***Revision date      yymmdd   (yymmdd)
C
C***Keywords
C***Author             Lengsfield, Byron (BRL)
C***Source             %W%   %G%
C
C***Purpose
C
C***Description
C
C***References
C
C***Routines called    (none)
C
C***End prologue
C
      implicit real*8(a-h,o-z)
      common /pdpack/ ipko, nhexo, ipkes
      data ipackl /'pack'/
c
c     positioning transformed integrals for reading
c     data set may or may not be in 'pack' format
c
c     rewind nfti
c     ipko = 1
c     read ( nfti ) ipack
c     if ( ipack .eq. ipackl ) go to 100
c     rewind nfti
 
c      call srew(nfti)
c
      ipko = 0
 100  return
      end
