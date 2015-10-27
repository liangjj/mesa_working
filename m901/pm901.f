*deck @(#)pm901.f	5.1  11/6/94
      subroutine pm901
c***begin prologue     pm901.f
c***date written       850601  
c***revision date      11/6/94      
c
c   3 december 1986   pws at lanl
c      changing 'namint' and iosys open to character.
c***keywords           
c***author             saxe,paul (lanl) 
c***source             @(#)pm901.f	5.1   11/6/94
c***purpose            
c***description
c     
c    
c
c***references
c
c***routines called
c
c***end prologue       pm901.f
c
c
      implicit integer (a-z)
c
      character*4096 ops
      character*2 mcorci
      character*8 citype
      character*128 gints
c
c     ----- record entries to distinguish from mcscf -----
      mcorci='ci'
      call iosys('write character mcorci to rwf',0,0,0,mcorci)
      citype='m901'
      call iosys('write character "ci used" to rwf',0,0,0,citype)
c
c     ----- recover the options string -----
c
      call iosys('read character options from rwf',-1,0,0,ops)
c
c     ----- call the ci routines -----
c
      call iosys('read character "guga integral filename" from rwf',
     $     0,0,0,gints)
      call iosys('open gints as old',0,0,0,gints)
      call mn901('guga integrals','gints')
c
c     ----- and exit with grace -----
c
      call iosys('close gints',0,0,0,' ')
      call chainx(0)
c
c
      return
      end
