*deck htree.f
c***begin prologue     htree
c***date written       961209   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           hartree, non-linear
c***author             schneider, barry (nsf)
c***source             diis
c***purpose            
c***                   
c***                   
c***references         
c
c***routines called    
c***end prologue       htree
      subroutine htree(ham0,ham,vec,pvec,vnl,eig,norm0,gamma,n)
      implicit integer (a-z)
      real*8 ham0, ham, vec, pvec, vnl, eig, norm0, gamma
      dimension ham0(n,n), vec(n), pvec(n,0:n-1), vnl(n), eig(n)
      dimension ham(n,n)
      common/io/inp, iout
      call copy(ham0,ham,n*n)
      call ebc(vnl,pvec,vec,n,n,1)
      do 10 i=1,n
         vnl(i) = vnl(i)*vnl(i)/(eig(i)*eig(i))
 10   continue
      call smul(vnl,vnl,gamma,n)
      do 20 i=1,n
         ham(i,i) = ham(i,i) + vnl(i)
 20   continue   
      return
      end       
