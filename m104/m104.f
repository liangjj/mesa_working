*deck @(#)m104.f	5.1  11/6/94
      program m104
c***begin prologue     m104.f
c***date written       yymmdd   (yymmdd)
c***revision date      11/6/94
c
c
c***keywords
c***author             lengsfield,byron (llnl)
c***source             @(#)pm104.f	5.1   11/6/94
c
c***purpose
c                      to reconfigure the rwf file for a kohn or r-matrix
c***                   calculation
c
c***description
c
c***references
c
c***routines called    (none)
c
c***end prologue       pm104.f
      implicit none
c     --- input variables -----
c     --- input arrays (unmodified) ---
c     --- input arrays (scratch) ---
c     --- output arrays ---
c     --- output variables ---
c     --- scratch arrays ---
c     --- local variables ---
      integer inp,iout
      integer nbf
      character*128 nmfile
      character*8 chrkey, scatyp, filtyp
      character*4096 ops
      real*8 rep
c
      common /io/ inp,iout 
c
      call drum
c     --- announce we're here.
      write(iout,*)' m104: '
      call iosys('read character options from rwf',-1,0,0,ops)
      scatyp=chrkey(ops,'scattering','none',' ')
c
c     since scatyp is a character*8 variable we have to examine only
c     the significant characters in making our decision.  trailing characters
c     need to be discarded.
c
      if(scatyp(1:4).eq.'none') then
         write(iout,1)
         call chainx(0)
         stop
      endif
      if(scatyp(1:4).eq.'kohn') then
         call iosys('read character "kohn filename" from rwf',
     #               -1,0,0,nmfile)
         filtyp='kohn'
         filtyp=filtyp(1:4)
      else if(scatyp(1:8).eq.'r-matrix') then
         call iosys('read character "r-matrix filename" from rwf',
     #               -1,0,0,nmfile)
         filtyp='rmtrx'
         filtyp=filtyp(1:5)
      endif
      call iosys('open '//filtyp//' as old',0,0,0,nmfile)
      call iosys('read integer "number of basis functions" from '
     $            //filtyp,1,nbf,0,' ')
      call iosys('write integer "number of basis functions" to rwf',
     $           1,nbf,0,' ')
c
      call iosys('read real "nuclear repulsion energy" from '//filtyp,
     $           1,rep,0,' ')
      call iosys('write real "nuclear repulsion energy" to rwf',
     $           1,rep,0,' ')
 
c
c
      call iosys('close '//filtyp,0,0,0,' ')
      call chainx(0)
c
      stop
 1    format(/,1x,'This link is only run to re-rwite the rwf file for '
     #            'scattering calculations. Exit') 
      end
