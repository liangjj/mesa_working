*deck hupdat.f
c***begin prologue     hupdat
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
c***end prologue       hupdat
      subroutine hupdat(h,htmp,b,btmp,vec,hvec,rhs,n,nrhs,nbeg,nend,m)
      implicit integer (a-z)
      real*8 h, htmp, b, btmp, vec, hvec, rhs, sdot
      dimension h(m,m), htmp(m,m), b(m,*), btmp(m,*), vec(n,*)
      dimension hvec(n,*), rhs(n,nrhs)
      common/io/inp, iout
      do 10 i=1,nend
         do 20 j=nbeg,nend
            h(i,j) = sdot(n,vec(1,i),1,hvec(1,j),1)
            h(j,i) = sdot(n,vec(1,j),1,hvec(1,i),1)
 20      continue
         do 30 j=1,nrhs   
            b(i,j) = sdot(n,vec(1,i),1,rhs(1,j),1)
 30      continue   
 10   continue
      do 40 i=1,nend
         do 50 j=1,i
            htmp(i,j) = h(i,j)
            htmp(j,i) = h(j,i)
 50      continue
         do 60 j=1,nrhs
            btmp(i,j) = b(i,j)
 60      continue   
 40   continue  
      return
      end       


