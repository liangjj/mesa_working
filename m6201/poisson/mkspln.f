*deck mkspln.f
c***begin prologue     mkspln
c***date written       930922   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           mkspln, link m6201, spline fit
c***author             schneider, barry (nsf)
c***source             m6201
c***purpose            calculate spline fit of partial wave functions
c***references         none
c
c***routines called
c***end prologue       mkspln
      subroutine mkspln (psilm,break,c,sc,mval,lval,n,order,
     1                   nbreak,nspln,strng)
      implicit integer (a-z)
      real*8 psilm, break, c, sc
      character*(*) strng
      dimension psilm(n,*), break(nbreak+1), c(order,nbreak)
      dimension sc(nspln,2)
      common /io/inp,iout
      call copy(sc(1,1),sc(1,2),nspln)
      locpf=0
      do 10 m=0,mval
         nm=2
         if (m.eq.0) then
             nm=1
         endif
         do 20 count=1,nm
            do 30 l=m,lval
               locpf=locpf+1
               call copy(sc(1,2),sc(1,1),nspln)
               call splcof(n,psilm(1,locpf),order,nbreak,break,c,sc)
               call iosys ('write real '//strng//' to lamdat '//
     1                     'without rewinding',order*nbreak,c,0,' ')                            
   30       continue 
   20    continue
   10 continue
      return
      end

