*deck tr2ply.f
c***begin prologue     tr2ply
c***date written       960723   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           basis, transformation
c***author             schneider, barry (nsf)
c***source             
c***purpose            transformation from dvr to polynomial basis.
c***                   
c***references         
c
c***routines called    
c***end prologue       tr2ply
      subroutine tr2ply(vec,trvec,tr,n,m,grid)
      implicit integer (a-z)
      real*8 vec, trvec, tr
      dimension vec(n,m), trvec(n,m), tr(n,n)
      character*2 itoc
      common/io/inp, iout 
      call iosys('read real "transformation matrix for grid '
     1           //itoc(grid)//'" from ham',n*n,tr,0,' ')        
      call ebc(trvec,tr,vec,n,n,m)
c      call iosys('write integer "number of vectors for grid '
c     1           //itoc(grid)//'" to ham',1,m,0,' ')     
      call iosys('write real "vectors for grid '//itoc(grid)
     1           //'" to ham',n*m,trvec,0,' ')     
      return
      end       
