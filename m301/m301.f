*deck  @(#)m301.f	1.3 7/30/91
      program m301
c
      implicit integer (a-z)
c
      character*4096 ops
      character*128 namint
      integer a
      real*8 z
      pointer (p,z(1)), (p,a(1))
      logical logkey
c
      data bufsiz/4096/
c
      common /io/ inp,iout
c
      call drum
c     see if we are to compute integrals.
      call iosys('read character options from rwf',-1,0,0,ops)
      if(logkey(ops,'int=reuse',.false.,' ').or.
     $   logkey(ops,'int=reuse1',.false.,' ').or.
     $   logkey(ops,'noints',.false.,' ')) then
c        restore the rwf file.
c        get some buffer space.
c         call manmem(0,idum,idum,'m301',0)
c         call manmem(bufsiz,p,ngot,'m301',0)
         call getmem(bufsiz,p,ngot,'m301',0)
c         call getscm(bufsiz,z(1),ngot,'m301: buffer',0)
         call iosys('read character "integral filename" from rwf',
     $               0,0,0,namint)
         call iosys('open ints as old',0,0,0,namint)
         call iosys('copy "packing index vector" from ints to rwf',
     $               bufsiz,a,0,'noerror')
         call iosys('copy "truncated number of basis functions"'
     $            //' from ints to rwf',bufsiz,a,0,'noerror')
         call iosys('close ints',0,0,0,' ')
c         call manmem(-ngot,p,idum,'m301',idum)
         call getmem(-ngot,p,idum,'m301',idum)
      else
c        call the routines to determine which functions to drop
         write(iout,*) (' m301:dropping functions')
         call mn301(ops)
      endif
c
c     ----- and exit with grace -----
c
      call chainx(0)
c
c
      stop
      end
