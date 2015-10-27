*deck conrmt
c***begin prologue     conrmt
c***date written       940213   (yymmdd)
c***revision date               (yymmdd)
c***keywords           r-matrix
c***author             schneider, barry (nsf)
c***source             
c***purpose            construct r-matrix
c***description        
c***references       
c
c***routines called
c***end prologue       conrmt
      subroutine conrmt(rmt,proj,eig,energy,nc,n,prnt)
      implicit integer (a-z)
      common /io/ inp, iout
      real*8 rmt, proj, eig, energy
      character*80 title
      logical prnt
      dimension rmt(nc,nc), proj(n,nc), eig(n)
      call rzero(rmt,nc*nc)
c      write(iout,*) eig
c      write(iout,*) proj
c      write(iout,*) energy
      do 10 i=1,nc
         do 20 j=1,i
            do 30 k=1,n
               rmt(i,j)=rmt(i,j)+proj(k,i)*proj(k,j)/(eig(k)-energy)
   30       continue
            rmt(i,j)=.5d0*rmt(i,j)
            rmt(j,i)=rmt(i,j)
   20    continue
   10 continue      
      if (prnt) then 
          title='r-matrix'
          call prntrm(title,rmt,nc,nc,nc,nc,iout)
      endif          
      return
      end
