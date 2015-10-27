*deck @(#)mysort.f	5.1  11/6/94
      subroutine mysort(flag,scr,iscr,ncor,len,label,labels,
     $     values,i1,c1,c2,c3)
c
c***begin prologue
c***date written       871022   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c
c***keywords
c***author             lengsfield, byron (brl)
c***source             @(#)mysort.f	5.1   11/6/94
c
c***purpose
c
c***description
c
c***references
c
c***routines called    (none)
c
c***end prologue
c
      implicit real*8 (a-h,o-z)
      dimension scr(*),iscr(*),label(*),labels(*),values(*)
c
      character*8 flag,c1,c2,file,unit,c3
c
      common /io/ inp,iout
c
c
      save lens,unit,file
c
      if(flag.eq.'start') then
c
         write(iout,*)' initialize sort space '
         write(iout,*)'  ncor  len ',ncor,len
c
         lens=len
         file=c1
         unit=c2
c
         if(len.gt.ncor) then
            call lnkerr(' out-of-core sort not implemented')
         endif
c
         call rzero(scr,len)
c
c
c
      elseif(flag.eq.'with bin') then
c
         write(iout,*)'  with bin  len = ',len
c
         write(iout,*)'  --- start scr(1) ',scr(1)
c
         do 1 i=1,len
            scr(label(i))=scr(label(i))+values(i)
 1       continue
c
         write(iout,*)'  --- end   scr(1) ',scr(1)
c
      else
c
         write(iout,*)' final sort '
         write(iout,*)' output is on file ',file
         write(iout,*)' file   is on unit ',unit
c
c      call printm(scr,5,1)
c
         call iosys('rewind mcscf_hessian on rwf',0,0,0,' ')
         call iosys('write real mcscf_hessian on rwf without rewinding',
     $        lens,scr,0,' ')
c
c
      endif
c
      return
      end
