      subroutine rdlocp (varr,vloc,ntchn,nptmx)
      implicit integer(a-z)
      dimension varr(ntchn,ntchn,nptmx)
      real *8 varr, fct
      character *8 vloc
      character *3 itoc
      data fct /2.d+00/
*
      nwords=ntchn*ntchn*nptmx
      call iosys ('open sc-mat as old',0,0,0,vloc)
      call iosys ('read real "lam locv-'//itoc(ntchn)//'" from '//
     1            'sc-mat',nwords,varr,0,0)
      call iosys ('rewind all on sc-mat read-and-write',0,0,0,0)
      call iosys ('close sc-mat',0,0,0,0)
      do 10 i=1,ntchn
      do 20 j=1,ntchn
      do 30 k=1,nptmx
      varr(i,j,k)=varr(i,j,k)*fct
   30 continue
   20 continue
   10 continue
      return
*
c
      end
