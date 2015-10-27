*deck setp2.f
c***begin prologue     setp2
c***date written       951229   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           two-dim
c***author             schneider, barry (nsf)
c***source             3-dim
c***purpose            form two dimensional index array.
c***                   
c***                   
c***references         
c
c***routines called    
c***end prologue       setp2
      subroutine setp2(ind,n,n1,n2,sym,prnt)
      implicit integer (a-z)
      logical prnt
      character*(*) sym
      dimension ind(n2,n1)
      common/io/inp, iout
      cnt=0 
      if(sym.eq.'unsymmetric') then
         do 10 i=1,n1
            do 20 j=1,n2
               cnt=cnt+1
               ind(j,i)=cnt
 20         continue
 10      continue   
         if(prnt) then
            write(iout,1)
            cnt=0
            do 30 i=1,n1
               do 40 j=1,n2
                  cnt=cnt+1 
                  write(iout,2) i, j, ind(j,i)
 40            continue   
 30         continue
         endif   
      elseif(sym.eq.'symmetric') then
         if(n1.ne.n2) then
            call lnkerr('error in symmetric case')
         endif
         do 50 i=1,n1
            do 60 j=1,i
               cnt=cnt+1 
               ind(i,j) = cnt
 60         continue
 50      continue   
         if(prnt) then
            write(iout,3)
            cnt=0
            do 70 i=1,n1
               do 80 j=1,i
                  cnt=cnt+1 
                  write(iout,2) i, j, ind(i,j)
 80            continue   
 70         continue
         endif   
      else
         call lnkerr('symmetry error')
      endif
      return
 1    format(/,1x,'hamiltonian indices for unsymmetric case')
 2    format('i = ',i3,1x,'j = ',i3,1x,'index = ',i3)
 3    format(/,1x,'hamiltonian indices for symmetric case')
      end       

