c \documentclass{article}
c \usepackage{graphicx}
c \setkeys{Gin}{width=\linewidth}
c \title{TRILIN: Read and Write Davidson Trial Vectors}
c \author{Barry I. Schneider}
c \date{}
c \def \<{\langle}
c \def \>{\rangle}
c \begin{document}
c \maketitle

*deck trilin.f
c***begin prologue     trilin
c***date written       010829   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           read in and write out trial vectors
c***author             schneider, barry (nsf)
c***source             
c***                   
c***references         
c
c***routines called    
c***end prologue       trilin
      subroutine trilin(vec,type,n,ntrial,maxvec,prnt)
      implicit integer (a-z)
      real*8 vec
      character*4 itoc
      character*80 title
      character*(*) type
      dimension vec(n,*)
      common/io/inp, iout
      call iosys('create real trials on rwf',n*ntrial,0,0,' ')
      start=0
      do while(start.lt.ntrial)
         end=min(start+maxvec,ntrial)
         start=start+1
         cnt=0
         do 10 i=start,end
            cnt=cnt+1
            call iosys('read real "'//type//itoc(i)//'" from rwf',
     1                  n,vec(1,cnt),0,' ')
 10      continue   
         call iosys('write real trials to rwf without rewinding',
     1               n*cnt,vec,0,' ')
         call iosys('rewind trials on rwf read-and-write',0,0,0,' ')
         if(prnt) then
            title='input vectors'
            call prntrm(title,vec,n,cnt,n,cnt,iout)
         endif
         start=end
      enddo
      return
      end

