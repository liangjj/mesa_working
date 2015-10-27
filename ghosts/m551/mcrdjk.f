*deck %W%  %G%
      subroutine mcrdjk(xjk,nblock,ndab,jkst,iocore)
C
C***Begin prologue     
C***Date written       871022   (yymmdd)
C***Revision date      yymmdd   (yymmdd)
C
C***Keywords           
C***Author             Lengsfield, Byron (BRL)
C***Source             %W%   %G%
C
C***Purpose            
C
C***Description        
C
C***References         
C
C***Routines called    (none)
C
C***End prologue       
C
c     implicit real*8(a-h,o-z)
cc
c
      dimension xjk(2)
      common / number / zero,pt5,one,two,four,eight
c---------------------------------------------------------------------c
c     this program reads blocks of coulomb(exchange) integrals
c     from a direct access data set
c     --each da block consists of nkl matrices and mrs integrals
c
c    xjk          coulomb(exchange) integrals
c    nrs          number of integrals in one coulomb(exchange) matrix
c    mrs          number of integrals in a da block
c    iijt         counter for the number matrices read
c    nkl          total number of matrices
c    mkl          number of matrices in a da block
c    jkst         starting da record
c    nblock       length of record to be read
c                 nblock = mkl * ( ncoljk * mrs )
c                 where ncoljk = (nrs-1)/mrs+1
c-----------------------------------------------------------------------
c
c.io      call wread(ndab,xjk,intowp(nblock),jkst,jkst)
c
      call iosys('read real mc_jk from mcscr',nblock,xjk,jkst,0)
      jkst=jkst+nblock
c
      return
      end
