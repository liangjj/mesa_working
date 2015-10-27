*deck cblslv.f
c***begin prologue     cblslv
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
c***end prologue       cblslv
      subroutine cblslv(vecin,vecout,temp,lufac,luind,n,nc,nvc)
      implicit integer (a-z)
      real*8 vecin, vecout
      complex*16 lufac, temp
      character*2 itoc
      dimension vecin(n,nc,2,nvc), vecout(n,nc,2,nvc)
      dimension lufac(*), luind(*), temp(n,nc,nvc)
      common/io/inp, iout
      call rv2cv(vecin,temp,n*nc,nvc)
      call iosys('read integer "number of LU blocks '//
     1           'per channel" from ham',1,ntrip,0,' ')
      do 10 ic=1,nc
         cnt=0
         do 20 trp=1,ntrip
            call iosys('read integer "block size for block-'//itoc(trp)
     1              //'channel-'//itoc(ic)//'" from ham',1,msize,0,' ')
            call iosys('read real "LU factors for block-'//itoc(trp)
     1              //'channel-'//itoc(ic)//'" from ham',
     2              2*msize*msize,lufac,0,' ')
            call iosys('read integer "LU array for block-'//itoc(trp)
     1                 //'channel-'//itoc(ic)//'" from ham',
     2                 msize,luind,0,' ') 
            do 30 i=1,nvc
               call cgetrs('n',msize,1,lufac,msize,
     1                      luind,temp(cnt+1,ic,i),n,info)
c               call cgesl(lufac,msize,msize,luind,temp(cnt+1,ic,i),0)
 30         continue   
            cnt=cnt+msize
 20      continue
 10   continue
      call cv2rv(vecout,temp,n*nc,nvc)
      return
      end       


