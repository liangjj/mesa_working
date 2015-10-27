*deck outvec.f
c***begin prologue     outvec
c***date written       980420   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           write out vectors
c***author             schneider, barry (nsf)
c***source             
c***                   
c***references         
c
c***routines called    
c***end prologue       outvec
      subroutine outvec(vec,type,n,nvc,prnt)
      implicit integer (a-z)
      real*8 vec
      character*4 itoc
      character*80 title
      character*(*) type
      dimension vec(n,nvc)
      common/io/inp, iout
      do 10 i=1,nvc
         call iosys('write real "'//type//itoc(i)//'" to rwf',
     1               n,vec(1,i),0,' ')
 10   continue
      if(prnt) then
         title='output vectors'
         call prntrm(title,vec,n,nvc,n,nvc,iout)
      endif
      return
      end

