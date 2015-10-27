*deck matprnt.f
c***begin prologue     matprnt
c***date written       951229   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           print routine
c***author             schneider, barry (nsf)
c***source
c***purpose            
c***                   
c***                   
c***description      
c***                 
c***                                                                       
c***                                                          
c***references         
c
c***routines called    
c***end prologue       matprnt

      subroutine matprnt (title,label,headr,a,nrow,ncol,ia,ib,ncmx)
      implicit integer (a-z)
      real*8 a, headr
      character*(*) label, title
      dimension a(ia,ib), headr(*)
      common/io/inp, iout 
      ntrip=ncol/ncmx
      final=ncmx
      left=ncol-ncmx*ntrip
      if(left.ne.0) then
         ntrip=ntrip+1
         final=left
      endif
      start=0
      do 10 i=1,ntrip-1
         ibeg=start+1
         iend=start+ncmx
         write(iout,1) label, ( headr(j),j=ibeg,iend)
         call prntfm(title,a(1,ibeg),nrow,ncmx,ia,ib,iout)
         start=start+ncmx
 10   continue   
      ibeg=start+1
      iend=start+final
      write(iout,1) label, ( headr(i),i=ibeg,iend)
      call prntfm(title,a(1,ibeg),nrow,final,ia,ib,iout)
      return
 1    format(/,a5,5e15.8)
      end       
