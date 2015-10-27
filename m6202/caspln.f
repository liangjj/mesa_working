*deck caspln.f
c***begin prologue     caspln
c***date written       930922   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           caspln, link m6201, spline fit
c***author             schneider, barry (nsf)
c***source             m6201
c***purpose            calculate partial wave functions from spline
c***references         none
c
c***routines called
c***end prologue       caspln
      subroutine caspln (psilm,r,break,c,ind,mval,lval,n,order,
     1                   nbreak,nspln,strng,prnt)
      implicit integer (a-z)
      real*8 psilm, break, c, r
      character*(*) strng
      logical prnt
      character*80 title
      dimension psilm(n,*), break(nbreak+1), c(order,nbreak)
      dimension ind(n), r(n)
      common /io/inp,iout
      call iosys ('rewind '//strng//' on lamdat read-and-write',
     1             0,0,0,' ')
      call fndbrk(r,break,ind,n,nbreak)
      locpf=0
      do 10 m=0,mval
         nm=2
         if (m.eq.0) then
             nm=1
         endif
         do 20 count=1,nm
            do 30 l=m,lval
               locpf=locpf+1
               call iosys ('read real '//strng//' from lamdat '//
     1                     'without rewinding',order*nbreak,c,0,' ')               
               call ppval(r,psilm(1,locpf),break,c,n,nbreak,order,ind,0)
   30       continue 
   20    continue
   10 continue
      if (prnt) then 
          title='fitted partial wave radial functions'
          call prntfm(title,psilm,n,locpf,n,locpf,iout) 
      endif           
      return
      end

