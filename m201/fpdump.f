*deck @(#)fpdump.f	5.1  11/6/94
      subroutine fpdump(nvar,pool0,pool1,delvar,yold,d1var,d2var,
     $                  d1vold,xi,h,vname)
c***begin prologue     fpdump.f
c***date written       yymmdd  
c***revision date      11/6/94      
c
c***keywords           
c***author             binkley, et al. (g82)
c***source             @(#)fpdump.f	5.1   11/6/94
c***purpose            
c***description
c     
c    
c
c***references
c
c***routines called
c
c***end prologue       fpdump.f
      implicit none
c     --- input variables -----
      integer nvar
c     --- input arrays (unmodified) ---
      character*(*) vname(nvar)
      real*8 pool0(nvar),pool1(nvar),delvar(nvar),yold(nvar)
      real*8 d1var(nvar),d2var(nvar),d1vold(nvar),xi(nvar)
      real*8 h(nvar,nvar)
c     --- input arrays (scratch) ---
c     --- output arrays ---
c     --- output variables ---
c     --- scratch arrays ---
c     --- local variables ---
      integer i,inp,iout
      integer fpcycl,mxcycl,mode,nstep,curvar,lambda
      real*8 convf,alpha,fzero,f1
      logical newcyc,inscal
c
      common/fpinf/ fpcycl,mxcycl,convf,mode,nstep,curvar,lambda,newcyc,
     $              inscal,alpha,fzero,f1(4)
      common/io/inp,iout
c
 1000 format(' fpdump, cycle ',i2,'.')
 1010 format(4x,'mxcycl    mode   nstep  curvar  lambda  newcyc',
     $          'inscal    alpha       fzero                f1')
 1020 format(2x,5(2x,i4),2(2x,l4),6e12.4)
 1030 format(8x,'variable',4x,'  pool0     ','  pool1     ',
     $                        '  delvar    ','  yold      ',
     $                        '  d1var     ','  d2var     ',
     $                        '  d1vold    ','  xi        ')
 1040 format(4x,a16,8d12.4)
 1050 format(2x,'h-matrix:')
c
c     --- dump the fletcher-powell information
      write(iout,1000)
      write(iout,1010)
      write(iout,1020) mxcycl,mode,nstep,curvar,lambda,newcyc,inscal,
     $                 alpha,fzero,(f1(i),i=1,4)
      write(iout,1030)
      do 10 i=1,nvar
         write(iout,1040) vname(i),pool0(i),pool1(i),delvar(i),yold(i),
     $                    d1var(i),d2var(i),d1vold(i),xi(i)
   10 continue
      write(iout,1050)
      call matprt(h,nvar,nvar,nvar,nvar,1,1,vname,vname,1,h,.false.)
c
c
      return
      end
