      subroutine setwpt (pt,wt,wts,ncst,nsts,nptmx,msize)
      implicit integer(a-z)
      real *8pt, wt, wts
      dimension pt(nptmx), wt(nptmx), wts(msize)
      dimension ncst(nsts)
      call iosys ('read real points from rwf',nptmx,pt,0,0)
      call iosys ('read real weights from rwf',nptmx,wt,0,0)
      icnt=0
      do 10 i=1,nsts
      nsch=ncst(i)
      do 10 j=1,nsch
      do 10 k=1,nptmx
      icnt=icnt+1
      wts(icnt)=wt(k)
   10 continue
      return
      end
