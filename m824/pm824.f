*deck @(#)pm824.f	5.1  11/6/94
      subroutine pm824(rcore,icore)
c
      implicit integer (a-z)
c
      character*4096 ops
      character*128 rdints,dints
      integer icore(*)
      real*8 rcore(*)
      logical mcscf,logkey
      dimension last(300)
      common /io/ inp,iout
c
      write(iout,11)
   11 format(1x,'m824: derivative lagrangians ')
c
c     ----- recover the options string -----
c
      call iosys('read character options from rwf',-1,0,0,ops)
c
c     ----- call the routines to make the supermatrices -----
c
      mcscf=logkey(ops,'mcscf',.false.,' ')
c
      call iosys('read integer "number of basis functions" from rwf',
     S           1,nbf,0,' ')
      nnp=nbf*(nbf+1)/2
c
      call iosys('read character "raw derivative integral filename"'
     $            //' from rwf',0,0,0,rdints)
      call iosys('open rdints as old',0,0,0,rdints)
c
      len=50*nnp
      call iosys('read character "derivative integral filename"'
     $            //' from rwf',0,0,0,dints)
      call iosys('open dints as unknown',len,0,0,dints)
c
      call iosys('read integer nder from rdints',1,nder,0,' ')
      if(nder.gt.300) call lnkerr(' m824: nder to big ')
      call iosys('read integer last from rdints',nder,last,0,' ')
c
      do 1 i=1,nder
         call mn333(rcore,icore,maxcor,i,last(i))
         if(mcscf) then
            call mn814(rcore,icore,i)
         else
            call mn1010(rcore,icore,i)
         end if
         call iosys('destroy scr',0,0,0,' ')
 1    continue
c
c
c     ----- and exit with grace -----
c
      call iosys('close rdints',0,0,0,' ')
      call iosys('close dints',0,0,0,' ')
      call chainx(0)
c
c
      stop
      end
