*deck scr2ham.f
c***begin prologue     scr2ham
c***date written       000921   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           hamiltonian
c***author             schneider, barry (nsf)
c***source             3-dim
c***purpose            move some buffered files.
c***                   
c***                   
c***references         
c
c***routines called    
c***end prologue       scr2ham
      subroutine scr2ham(ibuf,hbuf,lenbuf)
c
      implicit integer (a-z)
c
      real*8 hbuf
      character*40 bufr
      character*4 itoc
      dimension ibuf(2,lenbuf), hbuf(lenbuf)
      dimension bufr(2)
      common /io/ inp,iout
      bufr(1)='"contracted pq space buffers"'
      bufr(2)='"number of contracted pq space elements"'
      call iosys('read integer '//bufr(2)//' from scratch',
     1            1,nel,0,' ')
      call iosys('write integer '//bufr(2)//' to hamiltonian',
     1            1,nel,0,' ')
      call iosys('create integer '//bufr(1)//' on hamiltonian',
     1           -1,0,0,' ')
      trips=nel/lenbuf
      left=nel-trips*lenbuf
      do 20 nt=1,trips
         call iosys('read integer '//bufr(1)//' from scratch '//
     1              'without rewinding',2*lenbuf,ibuf,0,' ')
         call iosys('write integer '//bufr(1)//' to hamiltonian '//
     1              'without rewinding',2*lenbuf,ibuf,0,' ')
         call iosys('read integer '//bufr(1)//' from scratch '//
     1              'without rewinding',wptoin(lenbuf),hbuf,0,' ')
         call iosys('write integer '//bufr(1)//' to hamiltonian '//
     1              'without rewinding',wptoin(lenbuf),hbuf,0,' ')
 20   continue   
      if(left.ne.0) then
         call iosys('read integer '//bufr(1)//' from scratch '//
     1              'without rewinding',2*left,ibuf,0,' ')
         call iosys('write integer '//bufr(1)//' to hamiltonian '//
     1              'without rewinding',2*left,ibuf,0,' ')
         call iosys('read integer '//bufr(1)//' from scratch '//
     1              'without rewinding',wptoin(left),hbuf,0,' ')
         call iosys('write integer '//bufr(1)//' to hamiltonian '//
     1              'without rewinding',wptoin(left),hbuf,0,' ')
      endif   
      call iosys('endfile '//bufr(1)//' on hamiltonian',0,0,0,' ')
      call iosys('write integer '//bufr(2)//' to hamiltonian',
     1            1,nel,0,' ')
      call iosys('rewind all on hamiltonian read-and-write',0,0,0,' ')
      call iosys('destroy scratch',0,0,0,' ')
      write(iout,1) nel
      return
 1    format(/,1x,'number of non-zero contracted pq-space '
     1            'matrix elements = ',i5 )
      end
