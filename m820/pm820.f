*deck @@(#)m820.f	5.1  11/6/94
      program m820
c***begin prologue     m820
c***date written       871027   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c
c***keywords           guga integral ordering
c***author             saxe, paul (lanl)
c***source             @@(#)pm820.f	5.1   11/6/94
c
c***purpose            to sort mo ordered mo integrals to guga order.
c
c***description
c
c***references
c
c***routines called    (none)
c
c***end prologue       m820
c
      implicit integer (a-z)
c
      character*4096 ops
      character*128 tints,gints
      logical logkey
      integer icore(*)
      real*8 rcore(*)
c
      call drum
c
c     ----- recover the options string -----
c
      call iosys('read character options from rwf',-1,0,0,ops)
c
c     ----- call the routines to make the supermatrices -----
c           the integrals are read from tints; written to gints.
c
      call iosys('read character "transformed integral filename"'
     $         //' from rwf',0,0,0,tints)
      call iosys('open tints as old',0,0,0,tints)
c
      call iosys('read character "guga integral filename" from rwf',
     $            0,0,0,gints)
      len=0 
      if(logkey('ops','unit=ssd=gint',.false.,' ')) then
         call iosys('open gints as new on ssd',len,0,0,gints)
      else 
         call iosys('open gints as new',len,0,0,gints)
      endif
c
c     if we are sorting the tints file, we don't want to reorder
c     the orbitals in mn820. hence note argument 7.
      call mn820(rcore,icore,maxcor,'tints','gints',' ','ci')
c
c
c     ----- and exit with grace -----
c
      call iosys('close tints',0,0,0,' ')
      call iosys('close gints',0,0,0,' ')
      call chainx(0)
c
c
      stop
      end
