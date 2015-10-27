*deck @(#)spltst.f	1.1 9/8/91
c***begin prologue     m7000
c***date written       920525   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           m7000, link 7000, spline
c***author             schneider, b. (nsf)
c***source             m7004
c***purpose            test spline fitting
c***
c***references       
c
c***routines called    iosys, util and mdutil
c***end prologue       m7000
      subroutine spltst(f,df,g,dg,all,x,fs,xs,c,ind,n,nspln,
     1                  filnm,nfiles)
      implicit integer (a-z)
      real*8 f, df, g, dg, x, fs, xs, c, all, dum
      character*1 itoc
      character*4 rowlab, collab(6)
      character *(*) filnm
      dimension f(n), df(n), g(n), dg(n), x(n), all(n,5)
      dimension c(nspln), ind(nspln), fs(nspln), xs(nspln)
      common /io/ inp, iout
      collab(1)='r'
      collab(2)='f'
      collab(3)='df'
      collab(4)='g'
      collab(5)='dg'
      collab(6)='wron'
      call iosys ('rewind '//filnm//' on atomci read-and-write',0,
     1             0,0,' ')
c**********************************************************************c
c              find the array giving the index of the x point          c
c              in the original table which is to the left of the       c
c                           interpolating point.                       c
c**********************************************************************c
      call locate(xs,nspln,x,n,ind)
      do 200 l=1,nfiles
         call iosys ('read real '//filnm//' from atomci without '//
     1               'rewinding',nspln,fs,0,' ')
         call iosys ('read real '//filnm//' from atomci without '//
     1               'rewinding',nspln,c,0,' ')
         call splint(xs,fs,x,f,df,c,ind,nspln,n)
         call iosys ('read real '//filnm//' from atomci without '//
     1               'rewinding',nspln,fs,0,' ')
         call iosys ('read real '//filnm//' from atomci without '//
     1               'rewinding',nspln,c,0,' ')
         call splint(xs,fs,x,g,dg,c,ind,nspln,n)
         do 300 pt=1,n
            all(pt,6)=f(pt)*dg(pt)-g(pt)*df(pt)
  300    continue
         call matprt(all,n,6,n,6,0,1,rowlab,collab,0,dum,.false.)
  200 continue
      return
      end





















































