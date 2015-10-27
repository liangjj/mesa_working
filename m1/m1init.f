*deck @(#)m1init.f	5.1  11/6/94
      subroutine m1init(names,maxsiz)
c***begin prologue     m1init.f
c***date written       850601  yymmdd
c***revision date      11/6/94
c
c    1 december 1986   pws at lanl
c        changing 'names' to character and adding 'int' to the
c        replacement list
c
c***keywords           m1, i/o, files
c***author             martin, richard (lanl)
c***source             @(#)m1init.f	5.1   11/6/94
c***purpose            starts up link1.
c***description
c     call m1init(names,maxsiz)
c       names   list of file names.
c       maxsiz  maximum size defined by user.
c     m1init opens the input and output files and communicates
c     with the outside world regarding file name replacement and
c     memory management.
c***references
c***routines called    lxopen(mdutil)
c***end prologue       m1init.f
      implicit none
c     --- input variables -----
      integer maxsiz
c     --- input arrays (unmodified) ---
      character*(*) names(*)
c     --- input arrays (scratch) ---
c     --- output arrays ---
c     --- output variables ---
c     --- scratch arrays ---
c     --- local variables ---
      integer i,mchsiz
      integer inp,iout
c
      common/io/inp,iout
c
 1000 format(' main files/memory:')
 1010 format(5x,'inp:   ',a80,/,
     #       5x,'out:   ',a80,/,
     #       5x,'chk:   ',a48,5x,'dat:    ',a48,/,
     #       5x,'rwf:   ',a48,5x,'rint:   ',a48,/,
     #       5x,'int:   ',a48,5x,'tint:   ',a48,/,
     #       5x,'gint:  ',a48,5x,'rdint:  ',a48,/,
     #       5x,'dint:  ',a48,5x,'zint:   ',a48,/,
     #       5x,'ham:   ',a48,5x,'moden:  ',a48,/,
     #       5x,'aoden: ',a48,5x,'saoden: ',a48,/,
     #       5x,'gden:  ',a48,5x,'fci:    ',a48,/)
 1020 format(' other working files:')
 1030 format(5x,'kohn:    '  ,a48,5x,'kohndt: ',a48,/,
     #       5x,'grid:    '  ,a48,5x,'orbs:   ',a48,/,
     #       5x,'vstat:   '  ,a48,5x,'ylms:   ',a48,/,
     #       5x,'bessel:  '  ,a48,5x,'knints: ',a48,/,
     #       5x,'tmat:    '  ,a48,5x,'blktmat:',a48,/,
     #       5x,'optint:  '  ,a48,5x,'atomci: ',a48,/,
     #       5x,'fedvr:   '  ,a48,5x,'bec:    ',a48,/,
     #       5x,'hconfig: '  ,a48,5x,'rmtrx:  ',a48,/,
     #       5x,'tdse:    '  ,a48)
 1040 format(5x,'machine size: ',6x,i10)
c
c     --- open the input and output files, and communicate with the
c         outside world about various run parameters 
c         (file names,memory,..).
      call lxopen(names,mchsiz,maxsiz)
c
c     --- print the file replacement, memory requests and capabilities.
      write(iout,1000)
      write(iout,1010) (names(i),i=1,18)
      write(iout,1040) mchsiz
      write(iout,1020)
      write(iout,1030) (names(i),i=20,36)
c
c
      return
      end
