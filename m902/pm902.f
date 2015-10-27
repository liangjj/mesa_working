*deck @(#)pm902.f	5.1 11/6/94
      subroutine pm902(rcore,icore)
c***begin prologue     pm902
c***date written       850606   (yymmdd)  
c***revision date      880111   (yymmdd)
c
c 11 january 1988      bhl at brl
c    added option to compute the density without running the ci
c    m902=density is the appropiate flag .. this should allow one
c    to run the gugaci (m901) and compute the density with this program
c
c***keywords           
c***author             saxe,paul (lanl) 
c***source             @(#)pm902.f	5.1   11/6/94
c***purpose            
c         to construct hamiltonian and/or 1e and 2e density
c         from drt info
c***description
c     
c    
c
c***references
c
c***routines called
c
c***end prologue       pm902
      implicit integer (a-z)
c
      character*4096 ops
      character*2 mcorci
      character*8 citype
      character*128 gints
      character*128 gden
      integer icore(*)
      real*8 rcore(*)
      logical logkey
      common/io/inp,iout
c
c     ----- open the read-write file -----
c
      mcorci='ci'
      call iosys('write character mcorci to rwf',0,0,0,mcorci)
      citype='m902'
      call iosys('write character "ci used" to rwf',0,0,0,citype)
c
c     ----- recover the options string -----
c
      call iosys('read character options from rwf',-1,0,0,ops)
c
c     ----- call the ci routines -----
c
      maxcor=1
c
c     ----- open  the integral and density matrix file ----
c
      call iosys ('read character "guga integral filename" '//
     $            'from rwf',-1,0,0,gints)
      call iosys('open gints as old',0,0,0,gints)
      call iosys ('read character "guga density filename" '//
     $            'from rwf',-1,0,0,gden)
      call iosys('open gden as new',0,0,0,gden)
c
c     ----- ci section -----
c
      mcroot=0
      if(.not.logkey(ops,'m902=density',.false.,' ')) then
         write(iout,1)
  1      format(1x,'m902: hamiltonian construction')
         call mn902(icore,rcore,maxcor,'ci',0,0,'gints','gden',
     #              'ci',mcroot,'ci')
      end if
c
c     ----- density matrix segment -----
c
      if (logkey(ops,'ci=two-particle-density',.false.,' ').or.
     $    logkey(ops,'m902=density',.false.,' ')) then
         write(iout,2)
  2      format(1x,'m902: density construction')
         call mn902(icore,rcore,maxcor,'density',0,0,'gints','gden',
     $              'ci',mcroot,'ci')
      end if
c
c     ----- and exit with grace -----
c
      call iosys('close gints',0,0,0,' ')
      call iosys('close gden',0,0,0,' ')
      call chainx(0)
c
c
      stop
      end
