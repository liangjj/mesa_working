*deck @(#)pm840.f	5.1  11/6/94
      subroutine pm840(rcore,icore)
c***begin prologue     m840
c***date written       871119   (yymmdd)
c***revision date      910618   (yymmdd)
c
c   18 june    1991    rlm at lanl
c      reads ao density matrices from unit aoden and writes the 
c      group ordered matrices to soaden.
c   10 january 1988    bhl at brl
c      mcscf defaulted to .false. for multi-reference ci run
c
c***keywords           gradient dnesity matrix ordering
c***author             saxe, paul (lanl)
c***source             @(#)pm840.f	5.1   11/6/94
c
c***purpose            to sort the two-particle ao density matrix to
c  the order needed by the gradient codes.
c
c***description
c
c***references
c
c***routines called    (none)
c
c***end prologue       m840
c
      implicit integer (a-z)
c
      character*4096 ops
      character*128 aoden,saoden
      integer icore(*)
      real*8 rcore(*)
c
      maxcor=1
c
c     ----- recover the options string -----
c
      call iosys('read character options from rwf',-1,0,0,ops)
c
c     ----- call the routines to make the supermatrices -----
c
      call iosys('read character "ao density filename" from rwf',
     $            0,0,0,aoden)
      call iosys('open aoden as old',0,0,0,aoden)
c
      call iosys('read character "sao density filename" from rwf',
     $            0,0,0,saoden)
      call iosys('open saoden as new',0,0,0,saoden)
c
      call mn840(rcore,icore,maxcor)
c
      call iosys('close aoden',0,0,0,' ')
      call iosys('close saoden',0,0,0,' ')
c
c     ----- and exit with grace -----
c
      call chainx(0)
c
c
      stop
      end
