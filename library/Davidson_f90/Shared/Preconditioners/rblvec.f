*deck rblvec.f
c***begin prologue     rblvec
c***date written       960723   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           time development
c***author             schneider, barry (nsf)
c***source             
c***purpose            form new trial vectors using a block solve
c***                   strategy.
c***                                      
c***references         
c
c***routines called    
c***end prologue       rblvec
      subroutine rblvec(vecin,vecout,eig0,eig,u0,n,nvc)
      implicit integer (a-z)
      real*8 vecin, vecout, eig0, eig, u0
      character*2 itoc
      character*24 str
      dimension vecin(n,nvc), vecout(n,nvc), eig(nvc)
      dimension eig0(*), u0(*)
      common/io/inp, iout
      call iosys('read integer "number of blocks" from ham',1,
     1            ntrip,0,' ')
      cnt=0
      do 10 trp=1,ntrip
         str='block-'//itoc(trp)
         call iosys('read integer "block size for '//str//'" from '//
     1               'ham',1,msize,0,' ')
         call iosys('read real "eigenvalues for '//str//'" from '//
     1               'ham',msize,eig0,0,' ')
         call iosys('read real "transformation matrix for '//str
     1               //'" from ham',msize*msize,u0,0,' ') 
c
c        transform the input vectors to the representation diagonalizing 
c        the block.
c
         call ebtcxx(vecout(cnt+1,1),u0,vecin(cnt+1,1),msize,msize,nvc,
     1               n,msize,n)
c
c        energy scale the transformed vector
c
         call escale(vecout(cnt+1,1),vecin(cnt+1,1),eig,eig0,
     1               msize,nvc,n)
         call ebcxx(vecout(cnt+1,1),u0,vecin(cnt+1,1),msize,msize,nvc,
     1              n,msize,n)
         cnt=cnt+msize
 10   continue
      return
      end       


