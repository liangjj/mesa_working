      subroutine rhscal (gr1,gr2,rhs,ncst,jind,chnloc,nsts,nptmx,msize
     1 ,ntchn,ncmax,lplsmx,ops)
      implicit integer(a-z)
      logical logkey
      character *(*) ops
      real *8gr1, gr2, rhs
      dimension gr1(nptmx,ntchn), gr2(nptmx,ntchn), rhs(msize,ntchn)
      dimension ncst(nsts), jind(ncmax,nsts), chnloc(lplsmx,nsts)
      common /io/ inp, iout
      do 10 i=1,ntchn
      do 10 j=1,msize
   10 rhs(j,i)=0.d+00
      icnt=0
      do 20 is=1,nsts
      nsch=ncst(is)
      do 20 ls=1,nsch
      ncl=chnloc(jind(ls,is),is)
      do 20 ip=1,nptmx
      icnt=icnt+1
      rhs(icnt,ncl)=gr1(ip,ncl)*gr2(nptmx,ncl)
   20 continue
      if (logkey(ops,'print=lam=rhs',.false.,' ')) then
      write (iout,30)
      call matprt (rhs,msize,ntchn,msize,ntchn,0,0,0,0,0,0,0)
      endif
      return
c
   30 format (/,5x,'rhs')
      end
