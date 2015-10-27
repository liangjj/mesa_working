*deck rblslv.f
c***begin prologue     rblslv
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
c***end prologue       rblslv
      subroutine rblslv(vecin,vecout,lufac,luind,n,nc,nvc)
      implicit integer (a-z)
      real*8 vecin, vecout, lufac
      character*2 itoc
      dimension vecin(n,nc,nvc), vecout(n,nc,nvc), lufac(*), luind(*)
      common/io/inp, iout
      call copy(vecin,vecout,n*nc*nvc)
      call iosys('read integer "number of LU blocks '//
     1           'per channel" from ham',1,ntrip,0,' ')
      do 10 ic=1,nc
         cnt=0
         do 20 trp=1,ntrip
            call iosys('read integer "block size for block-'//itoc(trp)
     1              //'channel-'//itoc(ic)//'" from ham',1,msize,0,' ')
            call iosys('read real "LU factors for block-'//itoc(trp)
     1              //'channel-'//itoc(ic)//'" from ham',
     2              msize*msize,lufac,0,' ')
            call iosys('read integer "LU array for block-'//itoc(trp)
     1                 //'channel-'//itoc(ic)//'" from ham',
     2                 msize,luind,0,' ') 
            do 30 i=1,nvc
               call sgesl(lufac,msize,msize,luind,vecout(cnt+1,ic,i),0)
 30         continue   
            cnt=cnt+msize
 20      continue   
 10   continue
      return
      end       


