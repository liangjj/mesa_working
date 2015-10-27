      subroutine rdopt (vopt,optfle,rd,nopt,kpdiag,kepvnl,nsts,nfopt,ien
     1 ,ops)
      implicit integer(a-z)
      logical logkey
      dimension kpdiag(nsts,nsts), kepvnl(nfopt,nsts)
      dimension vopt(nfopt,nfopt), nopt(nsts)
      dimension rd(nfopt)
      common /io/ inp, iout
      character *80 rddag,rdodag
      character *3 itoc, ans
      character *8 optfle
      character *(*) ops
      real *8vopt, rd
      rddag(n1,n2)='read real "optdg'//itoc(n1)//itoc(n2)//'" from optic
     1al without rewinding'
      rdodag(n1,n2,n3)='read real "optodg'//itoc(n1)//itoc(n2)//itoc(n3)
     1 //'" from optical without rewinding'
      icpl=0
      if (logkey(ops,'zero-nlocv',.false.,' ')) icpl=1
      call iosys ('open optical as old',0,0,0,optfle)
      do 10 i=1,nfopt
      do 10 j=1,nfopt
   10 vopt(i,j)=0.d+00
      iread=ien
      call iosys ('does "optdg'//itoc(ien)//itoc(1)//'" exist on optical
     1 ',0,0,0,ans)
      if (ans.eq.'no') iread=1
      do 40 i=1,nsts
      norb=nopt(i)
      if (kpdiag(i,i).eq.0) then
      do 30 j=1,norb
      call iosys (rddag(iread,i),norb,rd,0,0)
      indi=kepvnl(j,i)
      do 20 k=1,j
      indj=kepvnl(k,i)
      vopt(indi,indj)=rd(k)
      vopt(indj,indi)=vopt(indi,indj)
   20 continue
   30 continue
      endif
   40 continue
      do 80 i=1,nsts
      i1=nopt(i)
      do 70 j=1,i
      if (i.eq.j) go to 70
      if (icpl.ne.0) go to 70
      if (kpdiag(i,j).eq.0) then
      j1=nopt(j)
      do 60 k=1,j1
      call iosys (rdodag(iread,i,j),i1,rd,0,0)
      indi=kepvnl(k,j)
      do 50 l=1,i1
      indj=kepvnl(l,i)
      vopt(indj,indi)=rd(l)
   50 vopt(indi,indj)=rd(l)
   60 continue
      endif
   70 continue
   80 continue
      call iosys ('rewind all on optical read-and-write',0,0,0,0)
      call iosys ('close optical',0,0,0,0)
      if (logkey(ops,'print=lam=optical-pot',.false.,' ')) then
      write (iout,90)
      call matprt (vopt,nfopt,nfopt,nfopt,nfopt,0,0,0,0,0,0,0)
      endif
      return
c
   90 format (/,5x,'optical potential')
      end
