*deck rmtrx
c***begin prologue     rmtrx
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
c***end prologue       rmtrx
      subroutine rmtrx(rmat,proj,eig,energy,nc,n,bypass,prnt)
      implicit integer (a-z)
      common /io/ inp, iout
      real*8 rmat, proj, eig, energy
      character*80 title
      logical prnt, bypass
      dimension rmat(nc,nc), proj(n,nc), eig(n)
      call rzero(rmat,nc*nc)
c      write(iout,*) eig
c      write(iout,*) proj
c      write(iout,*) energy
      if (.not.bypass) then
          do 10 i=1,nc
             do 20 j=1,i
                do 30 k=1,n
                   rmat(i,j)=rmat(i,j)+
     1                       proj(k,i)*proj(k,j)/(eig(k)-energy)
   30           continue
                rmat(i,j)=.5d0*rmat(i,j)
                rmat(j,i)=rmat(i,j)
   20        continue
   10     continue
      endif         
      if (prnt) then 
          title='r-matrix'
          call prntrm(title,rmat,nc,nc,nc,nc,iout)
      endif
      return
      end
