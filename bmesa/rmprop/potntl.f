*deck potntl
c***begin prologue     potntl
c***date written       940213   (yymmdd)
c***revision date               (yymmdd)
c***keywords           potential, matrix
c***author             schneider, barry (nsf)
c***source             
c***purpose            potential matrix elements
c***description        
c***references       
c
c***routines called
c***end prologue       potntl
      subroutine potntl(pt,vij,ec,vcij,range,type,nc,ntri)
      implicit integer (a-z)
      common /io/ inp, iout
      real*8 pt, vij, vcij, ec, range
      character*(*) type
      dimension vij(nc,nc), vcij(ntri), ec(nc)
      index=0
      do 10 i=1,nc
         do 20 j=1,i
            index=index+1
            call filpot(vij(i,j),vcij(index),pt,range,i,j,type)
            vij(j,i)=vij(i,j)
   20    continue
         vij(i,i)=vij(i,i)+ec(i)
   10 continue
      return
      end



