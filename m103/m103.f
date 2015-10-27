*deck @(#)pm103.f	5.1  11/6/94
      program m103
c***begin prologue     pm103.f
c***date written       850601   (yymmdd)
c***revision date      11/6/94
c
c        24 february 1989 by bhl at llnl
c           option to read maxsiz from ops
c***keywords           m103, link 1, $route, input
c***author             lengsfield, byron (llnl)
c***source             @(#)pm103.f	5.1   11/6/94
c***purpose            reset route options
c
c***description
c
c***references
c
c***routines called
c***end prologue       pm103.f
      implicit none
c     --- input variables -----
c     --- input arrays (unmodified) ---
c     --- input arrays (scratch) ---
c     --- output arrays ---
c     --- output variables ---
c     --- scratch arrays ---
c     --- local variables ---
      integer mxlnk,mxcard,mxlops
      parameter (mxlnk=200,mxcard=50)
      parameter (mxlops=4096)
c
      integer inp,iout
      integer lenops,maxsiz,i
      integer intkey
      character*4096 ops
      character*80 card
      logical dollar
c
      common/io/inp,iout
c
c
 1010 format(' title:')
 1020 format(5x,80a1)
 1030 format(' new options in force:')
c
c     --- retrieve the options string.
      if(dollar('$route',ops,card,inp)) then
      else
         call lnkerr('no $route section found.')
      endif
c     --- capitalize the input strings which will be parsed.
      call locase(ops,ops)
      call pakstr(ops,lenops)
      if(lenops.gt.mxlops) call lnkerr('option string too long.')
      call iosys('write character options to rwf',mxlops,0,0,ops)
c
      maxsiz=1000000
      maxsiz=intkey(ops,'maxsiz',maxsiz,' ')
      call iosys('write integer maxsiz to rwf',1,maxsiz,0,' ')
c
c     --- print the options string.
      write(iout,1030)
      write(iout,1020) (ops(i:i),i=1,lenops)
c
c     --- exit (gracefully?).
      call chainx(0)
c
c
      stop
      end
