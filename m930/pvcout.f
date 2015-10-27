*deck @(#)pvcout.f	5.1  11/6/94
      subroutine pvcout(vc,nspin,npvec,iu,filtyp)
      real*8 vc(nspin,*)
      character*(*) filtyp
      common/io/ inp,iout
c
      do 1 i=1,npvec
         do 2 j=1,nspin
            vc(j,i)=0.d0
   2   continue
         vc(i,i)=1.d0
   1  continue
c
      ntot=npvec*nspin
      call iosys('write integer "number of p-space vectors"'//
     $           ' to '//filtyp,1,npvec,0,' ')
      call iosys('write integer "number of p-space walks" to '//filtyp,
     $           1,nspin,0,' ')
      call iosys('write real "test vectors" to '//filtyp,
     $           ntot,vc,0,' ')
c
c
      return
      end
