*deck h1pac.f
c***begin prologue     h1pac
c***date written       000710   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           one-dim
c***author             schneider, barry (nsf)
c***source             3-dim
c***purpose            non zero hamiltonian matrix elements and indices 
c***                   for one dimensional hamiltonian. 
c***                   
c***                   
c***references         
c
c***routines called    
c***end prologue       h1pac
      subroutine h1pac(ham,v,hbuf,buf,diag,n,len,nonzro,prn)
      implicit integer (a-z)
      real*8 ham, v, hbuf, diag, zero
      logical prn
      character*80 title
      dimension ham(n,n), v(n), hbuf(len), buf(len,2), diag(n)
      data zero / 0.d0 /
      common/io/inp, iout
c
c     pack all non-zero, non-diagonal elements
c
      nonzro=0
      do 10 i=1,n
         do 20 j=1,i-1
            if(abs(ham(i,j)).gt.zero) then 
               nonzro=nonzro+1
               buf(nonzro,1) = i
               buf(nonzro,2) = j               
               hbuf(nonzro)=ham(i,j)
            endif               
 20      continue
 10   continue
c
c     diagonal elements contain both the kinetic energy and potential
c     matrix elements.
c 
      do 30 i=1,n
         diag(i)=ham(i,i)
 30   continue
      return
      end       


