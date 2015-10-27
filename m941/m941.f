*deck @(#)m901.f	5.1  11/6/94
      program m901
c***begin prologue     m901.f
c***date written       850601  
c***revision date      11/6/94      
c
c   3 december 1986   pws at lanl
c      changing 'namint' and iosys open to character.
c***keywords           
c***author             saxe,paul (lanl) 
c***source             @(#)m901.f	5.1   11/6/94
c***purpose            
c***description
c     
c    
c
c***references
c
c***routines called
c
c***end prologue       m901.f
c
c
      implicit integer (a-z)
c
      character*4096 ops
      character*2 mcorci
      character*8 citype
      character*128 gints
c      integer icore(maxcor)
c      real*8 rcore(*)
c
      call drum
c     ----- recover the options string -----
      call iosys('read character options from rwf',-1,0,0,ops)
c
c     ----- record entries to distinguish from mcscf -----
      mcorci='ci'
      call iosys('write character mcorci to rwf',0,0,0,mcorci)
      citype='m901'
      call iosys('write character "ci used" to rwf',0,0,0,citype)
c
c     ----- call the ci routines -----
c
      call iosys('read character "guga integral filename" from rwf',
     $     0,0,0,gints)
      call iosys('open gints as old',0,0,0,gints)
      call mn941(ops,'ci',0,0,'guga integrals','gints')
c
c     ----- and exit with grace -----
c
      call iosys('close gints',0,0,0,' ')
      call chainx(0)
c
c
      return
      end
