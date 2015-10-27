*deck prpac1.f
c***begin prologue     prpac1
c***date written       000710   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           one-dim
c***author             schneider, barry (nsf)
c***source             3-dim
c***purpose            print non zero hamiltonian matrix 
c***                   elements and indices for one dimensional hamiltonian. 
c***                   
c***                   
c***references         
c
c***routines called    
c***end prologue       prpac1
      subroutine prpac1(hbuf,buf,diag,n,len,nonzro)
      implicit integer (a-z)
      real*8 hbuf, diag
      character*80 title
      dimension hbuf(len), buf(len,2), diag(n)
      data zero / 0.d0 /
      common/io/inp, iout
      write(iout,1) nonzro
      write(iout,2)
      do 10 i=1,nonzro
         write(iout,3) buf(i,1), buf(i,2), hbuf(i)
 10   continue   
      write(iout,4) (diag(i),i=1,n)   
      return
 1    format(/,1x,'number of non-zero matrix elements = ',i5)
 2    format(/,5x,'    i   ',4x,'   j   ',3x,'  matrix element  ')
 3    format(5x,i5,6x,i5,6x,e15.8)
 4    format(/,1x,'diagonal elements = ',(/,5e15.8))
      end       


