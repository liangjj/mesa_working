*deck mkerr.f
c***begin prologue     mkerr
c***date written       960723   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           
c***author             schneider, barry (nsf)
c***source             
c***purpose            
c***                   
c***description          
c***references         
c
c***routines called    
c***end prologue       mkerr
      subroutine mkerr(err,obj1,obj2,iter,nobj,prnt)
      implicit integer (a-z)
      real*8 err, obj1, obj2
      character*80 title
      character*3 itoc
      logical prnt
      dimension err(*), obj1(*), obj2(*)
      common/io/ inp, iout
      do 10 i=1,nobj
         err(i) = obj1(i) - obj2(i)  
 10   continue
c      title='first object'
c      call prntrm(title,obj1,nobj,1,nobj,1,iout)
c      title='second object'
c      call prntrm(title,obj2,nobj,1,nobj,1,iout)            
      if(prnt) then
         title='error vector iteration = '//itoc(iter)
         call prntrm(title,err,nobj,1,nobj,1,iout)
      endif
      return
      end       

