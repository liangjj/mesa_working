*deck gusums.f
c***begin prologue     gusums
c***date written       951229   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           
c***author             schneider, barry (nsf)
c***source             
c***purpose            R-matrix sums for non-symmetric wavefunction.
c***                   
c***references         
c
c***routines called    
c***end prologue       gusums
      subroutine gusums(vec,srfx,srfy,sumx,sumy,n,nx,ny,nc,prn)
      implicit integer (a-z)
      real*8 vec, srfx, srfy, sumx, sumy
      character*80 title
      character*2 itoc
      logical prn
      dimension vec(ny,nx,n), sumx(nc,n), sumy(nc,n)
      common/io/inp, iout
c
c
      do 10 q=1,n
         do 20 k=1,nx
            sumx(k,q) = vec(ny,k,q)*srfy
 20      continue
         do 30 k=1,ny
            sumy(k,q)= vec(k,nx,q)*srfx
 30      continue   
 10   continue   
c      write(iout,1) srfx, srfy
c 1    format(/,1x,'srfx = ',e15.8,1x,'srfy = ',e15.8)
c      do 100 q=1,n
c         title='vector = '//itoc(q)
c         call prntrm(title,vec(1,1,q),ny,nx,ny,nx,iout)
c 100  continue   
      if(prn) then
         title='x sums'
         call prntrm(title,sumx,nx,n,nc,n,iout)      
        title='y sums'
         call prntrm(title,sumy,ny,n,nc,n,iout)      
      endif
      return      
      end       


