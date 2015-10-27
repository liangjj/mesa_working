*deck filodg
c***begin prologue     filodg
c***date written       940213   (yymmdd)
c***revision date               (yymmdd)
c***keywords           hamiltonian, matrix
c***author             schneider, barry (nsf)
c***source             
c***purpose            fill diagonal and off diagonal channel blocks
c***                   of the hamiltonian matrix.
c***description        
c***                   
c
c***routines called
c***end prologue       filodg
      subroutine filodg(h12,h21,v,ci,cj,ni,nj,n)
      implicit integer (a-z)
      common /io/ inp, iout
      real*8 h12, h21, v
      character*80 title
      dimension h12(n,n), h21(n,n), v(ni,nj)
c      write(iout,*) 'channel = ',ci,' channel = ',cj
c      title='v'
c      call prntrm(title,v,ni,nj,ni,nj,iout)
c      title='h12'
c      call prntrm(title,h12,ni,nj,n,n,iout)
c      title='h21'
c      call prntrm(title,h21,nj,ni,n,n,iout)
      if (ci.eq.cj) then
          do 10 i=1,ni
             do 20 j=1,i
                h12(i,j)=h12(i,j)+v(i,j)
                h12(j,i)=h12(i,j)
   20        continue
   10     continue
      else          
          do 30 i=1,ni
             do 40 j=1,nj
                h12(i,j)=h12(i,j)+v(i,j)
                h21(j,i)=h12(i,j)
   40        continue
   30     continue
      endif                                               
      return
      end
