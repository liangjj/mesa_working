*deck packh.f 
c***begin prologue     packh
c***date written       960718   (yymmdd)
c***revision date               (yymmdd)
c***keywords           pack, hamiltonian
c***author             schneider, b. i.(nsf)
c***source             
c***purpose            calculate a packed form of the DVR hamiltonian
c***                   
c
c***references       
c
c***routines called    iosys, util and mdutil
c***end prologue       packh
      subroutine packh(ham,d,buf,ibuf,nonz,n)
c
      implicit integer (a-z)
      real*8 ham, d, buf, eps
      dimension ham(n,n), d(n), buf(*), ibuf(2,*)
      data eps / 1.d-15 /
      common/io/inp, iout      
      posble=n*(n+1)/2    
      do 10 i=1,n
         d(i)=ham(i,i)
 10   continue
      nonz=0
      do 20 i=1,n
         do 30 j=1,i-1
            if(abs(ham(i,j)).gt.eps) then
               nonz=nonz+1
               buf(nonz)=ham(i,j)
               ibuf(1,nonz)=i
               ibuf(2,nonz)=j
            endif
 30      continue
 20   continue   
      count=nonz+n
      write(iout,1) count, posble
 1    format(/,1x,'number of non-zero matrix elements          = ',i10,
     1       /,1x,'possible number of non-zero matrix elements = ',i10)
      return
      end
