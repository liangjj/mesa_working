*deck @(#)schmab.f	5.1  11/6/94
      subroutine schmab(s,veca,vecb,ovab,t1,nbf,numa,numb,prnt,chklin)
c***begin prologue     schmab
c***date written       910603  (yymmdd)
c***revision date              (yymmdd)
c
c***keywords           vectors, orthogonalize
c***author             schneider, barry (lanl)
c***source
c***purpose            schmidt orthogonalize vector set b
c***                   to vector set a. set a is assumed
c***                   to be orthonormal.
c***                
c 
c***description        
 
c***references
c***routines called    aeqbc (mylib)
c                  
c                  
c***end prologue       schmab
      implicit integer (a-z)
c
      real*8 s(nbf,nbf), veca(nbf,numa), vecb(nbf,numb), ovab(numa,numb)
      real*8 t1(nbf,*)
      logical prnt, chklin
      common /io/ inp,iout
c
c
c   ----- s times vecb into t1 -----
c
      call aeqbc(t1,1,nbf,s,1,nbf,vecb,1,nbf,nbf,nbf,numb)
c
c   ----- transpose veca times t1 into ovab -----
c 
      call aeqbc(ovab,1,numa,veca,nbf,1,t1,1,nbf,numa,nbf,numb)
c
c    ----- reexpress vecb in the original basis -----
c
      call aeqbc(t1,1,nbf,veca,1,nbf,ovab,1,numa,nbf,numa,numb)
      do 10 i=1,numb
         do 20 j=1,nbf
            vecb(j,i)=vecb(j,i)-t1(j,i)
   20    continue
   10 continue
c
c     ----- print the vectors -----
c
      if (prnt) then
         write(iout,*) ' schmidt orthonormalized primitive set'
         call matout(vecb,nbf,numb,nbf,numb,iout)
      endif
c
c     ----- check orthogonality -----
      if (chklin) then
          call aeqbc(t1,1,nbf,s,1,nbf,vecb,1,nbf,nbf,nbf,numb)
          call aeqbc(ovab,1,numa,veca,nbf,1,t1,1,nbf,numa,nbf,numb)
          write(iout,*) ' schmidt overlap matrix'
          call matout(ovab,numa,numb,numa,numb,iout)
      endif
      return
      end
