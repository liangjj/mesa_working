*deck %W%  %G%
      subroutine pm551(cr,icr)
c
c***begin prologue     m551
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
c***end prologue       m551
c
      implicit real*8 (a-h,o-z)
c
      character*4096 ops
      character*8 prtflg,prtsav,mcscf
      character*128 namint,moden
      character*2 mcorci
      logical logkey
c
      integer icr(2)
c
      real*8 cr(1)
c
      common /io/ inp,iout
c
c
c
c
      write (iout,1)
 1    format(1x,'m551:  second order mcscf by b. lengsfield, b. liu',
     $          ' and m. yoshimine')
c
      mcorci='mc'
      call iosys('write character mcorci to rwf',0,0,0,mcorci)
c
c     ----- recover the options string -----
c
      call iosys('read character options from rwf',-1,0,0,ops)
      mcscf='mcscf'
      if (logkey(ops,'mcscf=cas',.false.,' ')) mcscf='casscf'
      if (logkey(ops,'mcscf=hf',.false.,' ')) mcscf='hf'
      if (logkey(ops,'hf=quadratic',.false.,' ')) mcscf='hf'
      if (mcscf.ne.'mcscf'.and.mcscf.ne.'casscf'
     $    .and.mcscf.ne.'hf') then
         call lnkerr('m551 must have "mcscf", "casscf" or "hf" '//
     $        'as an option')
      end if
c
c     ----- temporarily turn off printing as much as possible -----
c
      call iosys('read character "print flag" from rwf',-1,0,0,prtsav)
      prtflg='minimum'
      call iosys('write character "print flag" to rwf',0,0,0,prtflg)
c
c     ----- run the mcscf -----
c
      call iosys('read character "integral filename" from rwf',
     $            -1,0,0,namint)
      call iosys('open ints as old',0,0,0,namint)
c
c     open and write to the mo density file.
      call iosys('read character "mo density filename" from rwf',
     $            -1,0,0,moden)
      call iosys('open moden as new',0,0,0,moden)
c
      call iosys('write character mcscf_type to rwf',0,0,0,mcscf)
c
      call getscm(0,icr,ncore,'pm551',0)
      call mcdriv(cr,icr,ncore,prtflg)
c
c     ----- clean up -----
      call iosys('destroy mcscr',0,0,0,' ')
      call iosys('close ints',0,0,0,' ')
      call iosys('write character "print flag" to rwf',0,0,0,prtsav)
c
c
      call chainx(0)
c
      return
      end
