*deck pm1992
      subroutine pm1992(z,a)
c***begin prologue     m2
c***date written       920324   (yymmdd)
c***revision date               (yymmdd)
c
c
c***keywords
c***author             martin, richard  (lanl)
c***source             %W% %G% 
c***purpose
c                      determines direct access disk transfer rates 
c***description
c***references
c
c***routines called
c
c***end prologue       m1992
c
      implicit integer (a-z)
c
      integer a(1)
      real*8 z(1)
c     
      real*8 t1,t2,time
      real*8 rate
c
      character*4 itoc
      parameter (nrec=1,length=1000000)
      common/io/inp,iout
c
c     write out a few records and time it.
      call rzero(z(1),length)
c**** need a reliable timing indicator here
      call secnds(t1)
      do 10 i=1,nrec
	 call iosys('write real dummy'//itoc(i)//' on rwf',
     $               length,z(1),0,' ')
   10 continue
      call secnds(t2)
      time=t2-t1
      write(iout,*) 'write test:  nrecords length total time',
     $               nrec,length,nrec*length,time
      total=nrec*length
      rate=total*wptbyt(1)/time
      write(iout,*) 'rate(mbytes/sec):',rate
c
c     read them back in and check that.
      call secnds(t1)
      do 20 i=1,nrec
	 call iosys('read real dummy'//itoc(i)//' from rwf',
     $	             length,z(1),0,' ')
   20 continue
      call secnds(t2)
      time=t2-t1
      rate=total*wptbyt(1)/time
      write(iout,*) 'read test'
      write(iout,*) rate
c
      call chainx(0)
c
c
      stop
      end
