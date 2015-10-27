*deck hinit.f
c***begin prologue     hinit
c***date written       960723   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           time development
c***author             schneider, barry (nsf)
c***source             
c***purpose            
c***                                      
c***references         
c
c***routines called    
c***end prologue       hinit
      subroutine hinit(h,htmp,b,btmp,vec,hvec,rhs,n,nrhs,nend,m)
      implicit integer (a-z)
      real*8 h, htmp, b, btmp, vec, hvec, rhs, sdot
      dimension h(m,m), htmp(m,m), b(m,*), btmp(m,*), vec(n,*)
      dimension hvec(n,*), rhs(n,nrhs)
      common/io/inp, iout
      do 10 i=1,nend
         do 20 j=1,i
            h(i,j) = sdot(n,vec(1,i),1,hvec(1,j),1)
            h(j,i) = sdot(n,hvec(1,i),1,vec(1,j),1)
            htmp(i,j) = h(i,j)
            htmp(j,i) = h(j,i)
 20      continue
         do 30 j=1,nrhs   
            b(i,j) = sdot(n,vec(1,i),1,rhs(1,j),1)
            btmp(i,j) = b(i,j)
 30      continue	    
 10   continue
      return
      end       


