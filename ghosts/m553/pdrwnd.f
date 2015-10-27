      subroutine pdrwnd ( nfti )
c
c***begin prologue
c***date written       871022   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c
c***keywords
c***author             lengsfield, byron (brl)
c***source             %W%   %G%
c
c***purpose
c
c***description
c
c***references
c
c***routines called    (none)
c
c***end prologue
c
c     implicit real*8(a-h,o-p,r-z),integer*2(q)
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
