*deck @(#)trans.f	5.1  11/6/94
      subroutine trans(half,conint,ic,id,nab,fc,fd,t,lent,c,d)
c
      implicit integer (a-z)
c
c     real half(lenblk,fa,fb,ic,id),conint(lenblk,fa,fb,fc,fd)
      real*8 half(nab,ic,id),conint(nab,fc,fd)
      real*8 t(lent),c(ic,fc),d(id,fd)
      common/io/inp,iout
c
c     ----- start timing -----
c
c
c      write(iout,*) 'zeroing conint'
c      write(iout,*) 'lent',lent
      call rzero(conint,nab*fc*fd)
c
      n=min(nab,lent/(fc*id))
      npass=(nab+n-1)/n
      maxn=0
      do 100 pass=1,npass
         minn=maxn+1
         if (minn.gt.nab) go to 101
         maxn=min(nab,maxn+n)
         n=maxn-minn+1
c
c         write(iout,*) 'zeroing t',n*fc*id
         call rzero(t,n*fc*id)
c         write(iout,*) 'finished zeroing t'
c
         do 10 fcpt=1,fc
            do 9 icpt=1,ic
c               write(iout,*) 'fcpt icpt',fcpt,icpt
               if (c(icpt,fcpt).ne.0.0d+00) then
                  do 8 idpt=1,id
                     fcid=n*(fcpt-1+fc*(idpt-1))+1
c                     write(iout,*) 'fcid',fcid
c                     write(iout,*) 'minn icpt idpt',minn,icpt,idpt
                     call saxpy(n,c(icpt,fcpt),half(minn,icpt,idpt),1,
     #                          t(fcid),1)
    8             continue
               end if
    9       continue
   10    continue
c       write(iout,*) 'finished 10 loop'
c
         do 20 fdpt=1,fd
            do 19 idpt=1,id
               if (d(idpt,fdpt).ne.0.0d+00) then
                  do 18 fcpt=1,fc
                     fcid=n*(fcpt-1+fc*(idpt-1))+1
                     call saxpy(n,d(idpt,fdpt),t(fcid),1,
     #                          conint(minn,fcpt,fdpt),1)
   18             continue
               end if
   19       continue
   20    continue
c       write(iout,*) 'finished 20 loop'
c
  100 continue
  101 continue
c
c     do 100 n=1,lenblk
c        do 90 i=1,fa
c           do 80 j=1,fb
c
c              kl=0
c              do 2 l=1,id
c                 do 1 k=1,ic
c                    kl=kl+1
c                    t1(kl)=half(n,i,j,k,l)
c   1             continue
c   2          continue
c
c              call ebtc(t2,c,t1,fc,ic,id)
c              call ebc(t1,t2,d,fc,id,fd)
c
c              kl=0
c              do 4 l=1,fd
c                 do 3 k=1,fc
c                    kl=kl+1
c                    conint(n,i,j,k,l)=t1(kl)
c   3             continue
c   4          continue
c
c  80       continue
c  90    continue
c 100 continue
c
c     ----- stop timing -----
c
c
c
      return
      end
