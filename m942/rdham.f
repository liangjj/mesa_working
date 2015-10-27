*deck rdham.f
c***begin prologue     rdham
c***date written       960615   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           hamiltonian
c***author             schneider, barry (nsf)
c***source             3-dim
c***purpose            print non-zero matrix elements
c***                   of hamiltonian and their indices.
c***                   
c***                   
c***references         
c
c***routines called    
c***end prologue       rdham
      subroutine rdham(ibuf,hbuf,header,lenbuf,nonzro)
      implicit integer (a-z)
      real*8 hbuf
      character*(*) header
      dimension hbuf(lenbuf), ibuf(2,lenbuf)
      common/io/inp, iout 
      if(nonzro.eq.0) then
         return
      endif
      trips = nonzro/lenbuf
      left = nonzro - trips*lenbuf
      call iosys('rewind '//header//' on hamiltonian '//
     #           'read-and-write',0,0,0,' ')
      do 10 i=1,trips
         call iosys('read integer '//header//' from hamiltonian '//
     #              'without rewinding',2*lenbuf,ibuf,0,' ')
         call iosys('read integer '//header//' from hamiltonian '//
     #              'without rewinding',wptoin(lenbuf),hbuf,0,' ')
 20      continue   
 10   continue   
      if(left.ne.0) then
         call iosys('read integer '//header//' from hamiltonian '//
     #              'without rewinding',2*left,ibuf,0,' ')
         call iosys('read integer '//header//' from hamiltonian '//
     #              'without rewinding', wptoin(left),hbuf,0,' ')
         do 30 j=1,left
            write(iout,*) ibuf(1,j), ibuf(2,j), hbuf(j)
 30      continue   
      endif   
      call iosys('rewind all on hamiltonian read-and-write',0,0,0,' ')
      return
      end       


















