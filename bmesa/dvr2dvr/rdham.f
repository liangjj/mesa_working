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
      subroutine rdham(hbuf,ibuf,diag,q3,lenbuf,matdim,grid)
      implicit integer (a-z)
      real*8 hbuf, diag
      character*2 itoc
      character*80 title
      dimension hbuf(lenbuf), ibuf(2,lenbuf), diag(matdim)
      common/io/inp, iout 
      call iosys('read integer "'//itoc(q3)//'d number of elements '
     1           //'for grid '//itoc(grid)//'" from ham',1,ntot,0,' ')
      call iosys('read real "'//itoc(q3)//'d diagonals '
     1          //'for grid '//itoc(grid)//'" from ham',matdim,
     2             diag,0,' ')
      title='diagonal hamiltonian matrix elements'
      call prntrm(title,diag,matdim,1,matdim,1,iout)
      trips=ntot/lenbuf
      left=ntot-trips*lenbuf
      write(iout,1)
      write(iout,2)
      do 20 i=1,trips
         call iosys('read integer "'//itoc(q3)//
     1              'd buffers for grid '//itoc(grid)
     2              //'" from ham without '//
     3                'rewinding',2*lenbuf,ibuf,0,' ')
         call iosys('read integer "'//itoc(q3)//
     1              'd buffers for grid '//itoc(grid)
     2              //'" from ham without '//
     3                'rewinding',wptoin(lenbuf),hbuf,0,' ')
         do 30 j=1,lenbuf
            write(iout,3) ibuf(1,j), ibuf(2,j), hbuf(j)
 30      continue   
 20   continue   
      if(left.ne.0) then
         call iosys('read integer "'//itoc(q3)//
     1              'd buffers for grid '//itoc(grid)
     2              //'" from ham without '//
     3                'rewinding',2*left,ibuf,0,' ')
         call iosys('read integer "'//itoc(q3)//
     1              'd buffers for grid '//itoc(grid)
     2              //'" from ham without '//
     3                'rewinding',wptoin(left),hbuf,0,' ')
         do 40 j=1,left
            write(iout,3) ibuf(1,j),ibuf(2,j), hbuf(j)
 40      continue   
      endif
      return
 1    format(/,1x,'non-zero off diagonal matrix elements')
 2    format(/,5x,'    i   ',4x,'   j   ',4x,'  matrix element  ')
 3    format(5x,i5,6x,i5,6x,e15.8)
      end       






