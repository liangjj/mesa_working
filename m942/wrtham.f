*deck wrtham.f
c***begin prologue     wrtham
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
c***end prologue       wrtham
      subroutine wrtham(ibuf,hbuf,lenbuf,nonzro)
      implicit integer (a-z)
      real*8 hbuf
      dimension hbuf(lenbuf), ibuf(2,lenbuf)
      common/io/inp, iout 
      do 10 i=1,nonzro
         write(iout,*) ibuf(1,i), ibuf(2,i), hbuf(i)
 10   continue   
      return
      end       


















