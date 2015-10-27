*deck @(#)ectb.f	5.1  11/6/94
      subroutine ectb(ec,nnp,nkwks,up,iup,wpti4,
     #                ntops,dn,idn,nbots,jdnnwk,kdnnwk,
     #                offset,c,nwks,asort,jsort,jnsort,aov,jov,
     #                kov,jnov,maxov,mink,maxk,nij,ijpt)
c
      implicit integer (a-z)
c
c     ----- external arrays modified -----
c
      real*8 ec(nnp,nkwks)
c
c     ----- external arrays not modified -----
c           note that up,iup  and dn,idn are equivalenced in the
c           calling routine.
c
      real*8 c(nwks),up(4,ntops),dn(4,nbots)
      integer iup(wpti4,ntops),idn(wpti4,nbots),ijpt(nij)
c
c     ----- external scratch arrays -----
c
      real*8 asort(kdnnwk),aov(maxov)
      integer jsort(kdnnwk),jnsort(kdnnwk),jov(maxov),kov(maxov)
      integer jnov(maxov)
c
c     ----- external scalars not modified -----
c
      integer nnp,nkwks,ntops,nbots,jdnnwk,kdnnwk,nwks,offset
c
c     ----- local scalars -----
c
      integer ii,ktop,jtop,ij,j,k
c
c     ----- do loop variables -----
c
      integer top,bot,upwk,dnwk
c
      common /io/ inp,iout
c
c     ----- sort the bottom segments -----
c
      call izero(jsort,kdnnwk)
      call rzero(asort,kdnnwk)
      do 18 i=1,kdnnwk
         jnsort(i)=1
   18 continue
c
      n=0
      do 20 bot=1,nbots
         j=idn(2,bot)
         k=idn(1,bot)
         do 19 dnwk=1,1
            j=j+1
            k=k+1
            if (jsort(k).ne.0) then
               n=n+1
               if (n.gt.maxov) call lnkerr('too many overflows')
               jov(n)=j
               kov(n)=k
               jnov(n)=idn(3,bot)
               aov(n)=dn(4,bot)
            else
               jsort(k)=j
               jnsort(k)=idn(3,bot)
               asort(k)=dn(4,bot)
            end if
   19    continue
   20 continue
c
c
c     ----- form e(ij,k,j) * c(j) -----
c
      do 10 top=1,ntops
         ii=iup(3,top)*(iup(3,top)-1)/2
         ktop=(iup(1,top)-mink+1)*kdnnwk
         jtop=iup(2,top)*jdnnwk+offset
         kupwk=iup(1,top)
c
         do 40 upwk=1,1
c
            kupwk=kupwk+1
            if (kupwk.gt.maxk) go to 10
            if (kupwk.lt.mink) go to 36
c
cdir$ ivdep
            do 30 k=1,kdnnwk
               ij=ijpt(ii+jnsort(k))
               ec(ij,k+ktop)=ec(ij,k+ktop)+
     #                         up(4,top)*asort(k)*c(jsort(k)+jtop)
   30       continue
c
            do 35 i=1,n
               ij=ijpt(ii+jnov(i))
               ec(ij,kov(i)+ktop)=ec(ij,kov(i)+ktop)+
     #                              up(4,top)*aov(i)*c(jov(i)+jtop)
   35       continue
c
   36       continue
            ktop=ktop+kdnnwk
            jtop=jtop+jdnnwk
   40    continue
c
   10 continue
c
c
      return
      end
