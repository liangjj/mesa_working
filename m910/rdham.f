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
      subroutine rdham(hbuf,ibuf,diag,lenbuf,n,title)
      implicit integer (a-z)
      real*8 hbuf, diag
      character*(*) title
      dimension hbuf(lenbuf), ibuf(2,lenbuf), diag(n)
      dimension title(3)
      common/io/inp, iout 
      if(n.gt.0) then
         call iosys('read real '//title(3)//' from hamiltonian',
     1               n,diag,0,' ')
         write(iout,1) (diag(i),i=1,n)
      endif
      call iosys('read integer '//title(2)//' from hamiltonian',
     1            1,nonzro,0,' ')
      if(nonzro.eq.0) then
         return
      endif
      trips = nonzro/lenbuf
      left = nonzro - trips*lenbuf
      write(iout,2)
      write(iout,3)
      do 10 i=1,trips
         call iosys('read integer '//title(1)//
     1              ' from hamiltonian without rewinding',
     2                2*lenbuf,ibuf,0,' ')
         call iosys('read integer '//title(1)//
     1              ' from hamiltonian without rewinding',
     2                wptoin(lenbuf),hbuf,0,' ')
         do 20 j=1,lenbuf
            write(iout,4) ibuf(1,j),ibuf(2,j), hbuf(j)
 20      continue   
 10   continue   
      if(left.ne.0) then
         call iosys('read integer '//title(1)//
     1              ' from hamiltonian without rewinding',
     2                2*left,ibuf,0,' ')
         call iosys('read integer '//title(1)//
     1              ' from hamiltonian without rewinding',
     2                wptoin(left),hbuf,0,' ')
         do 30 j=1,left
            write(iout,4) ibuf(1,j),ibuf(2,j), hbuf(j)
 30      continue
      endif   
      call iosys('rewind all on hamiltonian read-and-write',0,0,0,' ')
      return
 1    format(/,1x,'diagonal elements = ',(/,5e15.8))
 2    format(/,1x,'non-zero matrix elements')
 3    format(/,5x,'    i   ',4x,'   j   ',3x,'  matrix element  ')
 4    format(5x,i5,6x,i5,6x,e15.8)
      end       


















