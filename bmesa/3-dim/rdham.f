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
      subroutine rdham(hbuf,ibuf,diag,q3,lenbuf,matdim,title)
      implicit integer (a-z)
      real*8 hbuf, diag
      character*(*) title
      character*1 itoc
      dimension hbuf(lenbuf), ibuf(2,lenbuf), diag(matdim)
      common/io/inp, iout 
      call iosys('read integer "'//itoc(q3)//'d number of elements '
     1           //title//'" from ham',1,ntot,0,' ')
      call iosys('read real "'//itoc(q3)//'d diagonals '
     1          //title//'" from ham',matdim,diag,0,' ')
      write(iout,1)
      write(iout,2)
      do 10 i=1,matdim
         write(iout,3) i, diag(i)
 10   continue   
      trips=ntot/lenbuf
      left=ntot-trips*lenbuf
      write(iout,4)
      write(iout,5)
      do 20 i=1,trips
         call iosys('read integer "'//itoc(q3)//
     1              'd buffers '//title//'" from ham '//
     2              'without rewinding',2*lenbuf,
     3               ibuf,0,' ')
         call iosys('read integer "'//itoc(q3)//
     1              'd buffers '//title//'" from ham '//
     2              'without rewinding',wptoin(lenbuf),
     3               hbuf,0,' ')
         do 30 j=1,lenbuf
            write(iout,6) ibuf(1,j),ibuf(2,j), hbuf(j)
 30      continue   
 20   continue   
      if(left.ne.0) then
         call iosys('read integer "'//itoc(q3)//
     1              'd buffers '//title//'" from ham '//
     2              'without rewinding',2*left,
     3               ibuf,0,' ')
         call iosys('read integer "'//itoc(q3)//
     1              'd buffers '//title//'" from ham '//
     2              'without rewinding',wptoin(left),
     3               hbuf,0,' ')
         do 40 j=1,left
            write(iout,6) ibuf(1,j),ibuf(2,j), hbuf(j)
 40      continue   
      endif
      return
 1    format(/,1x,'diagonal hamiltonian matrix elements')
 2    format(/,5x,'  index ',3x,'  matrix element  ')
 3    format(5x,i5,6x,e15.8)
 4    format(/,1x,'non-zero off diagonal matrix elements')
 5    format(/,5x,'    i   ',4x,'   j   ',3x,'  matrix element  ')
 6    format(5x,i5,6x,i5,6x,e15.8)
      end       






