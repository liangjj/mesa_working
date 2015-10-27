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
      subroutine rdham(hbuf,ibuf,diag,lenbuf,nonzro,incore,n)
      implicit integer (a-z)
      real*8 hbuf, diag
      logical incore
      dimension hbuf(lenbuf), ibuf(lenbuf,2), diag(n)
      common/io/inp, iout 
      if(incore) then
         do 10 i=1,nonzro
            write(iout,3) ibuf(i,1),ibuf(i,2), hbuf(i)
 10      continue   
      else
         call iosys('read integer "number of non-zero '//
     1              'matrix elements" from ham',1,nonzro,0,' ')
         trips=nonzro/lenbuf
         left=nonzro-trips*lenbuf
         write(iout,1)
         write(iout,2)
         do 20 i=1,trips
            call iosys('read integer "hamiltonian buffers" from ham '//
     1                 'without rewinding',2*lenbuf,ibuf,0,' ')
            call iosys('read integer "hamiltonian buffers" from ham '//
     1                 'without rewinding',wptoin(lenbuf),hbuf,0,' ')
            do 30 j=1,lenbuf
               write(iout,3) ibuf(j,1),ibuf(j,2), hbuf(j)
 30         continue   
 20      continue   
         if(left.ne.0) then
            call iosys('read integer "hamiltonian buffers" from ham '//
     1                 'without rewinding',2*left,ibuf,0,' ')
            call iosys('read integer "hamiltonian buffers" from ham '//
     1                 'without rewinding',wptoin(left),hbuf,0,' ')
            do 40 j=1,left
               write(iout,3) ibuf(j,1),ibuf(j,2), hbuf(j)
 40         continue   
         endif
         call iosys('rewind all on ham read-and-write',0,0,0,' ')
      endif
      write(iout,4) (diag(i),i=1,n)
      return
 1    format(/,1x,'non-zero matrix elements')
 2    format(/,5x,'    i   ',4x,'   j   ',3x,'  matrix element  ')
 3    format(5x,i5,6x,i5,6x,e15.8)
 4    format(/,1x,'diagonal elements = ',(/,5e15.8))
      end       


















