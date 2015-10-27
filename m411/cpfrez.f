*deck @(#)cpfrez.f	5.1  11/6/94
      subroutine cpfrez(vecold,eigold,symold,vecnew,eignew,
     1                  symnew,index,nbf,freeze,sym)
c
c***begin prologue     cpfrez
c***date written       910606  (yymmdd)
c***revision date              (yymmdd)
c
c***keywords           vectors, transfer
c***author             schneider, barry (lanl)
c***source
c***purpose            to copy a subset of an old vector
c***                   set (vecold) to a new vector set (vecnew)
c***                   and to then return it to the initial set
c***                   in ordinal order
c 
c***description
c***references
c***routines called    scopy, icopy
c                      
c***end prologue       cpfrez
      implicit integer (a-z)
      real*8 vecold, eigold, vecnew, eignew
      logical sym
      dimension vecold(nbf,nbf), vecnew(nbf,freeze), eigold(nbf)
      dimension eignew(freeze), symold(nbf), symnew(freeze), index(nbf)
      do 10 nfreze=1,freeze  
         call scopy(nbf,vecold(1,index(nfreze)),1,vecnew(1,nfreze),1)
         eignew(nfreze)=eigold(index(nfreze))
   10 continue
      call scopy(freeze,eignew,1,eigold,1)
      if (sym) then
          do 20 nfreze=1,freeze
             symnew(nfreze)=symold(index(nfreze))
   20     continue 
          call icopy(symnew,symold,freeze)
      endif
      return
      end   
