*deck @(#)sort32.f	5.1  11/6/94
      subroutine sort32(labels,values,acore,lenbuf,num,nnp,ints,iout,
     #                 maxcor,vals,labs,bins,lenscr,ijv,klv,
     #                 iork,jorl,prnt,nder,last)
c
c   nb.  iork and jorl are equivalenced to labs by the calling sequence
c
      implicit integer (a-z)
c
      logical prnt
      logical debug
      real*8 values(lenbuf)
      integer acore(maxcor)
      character*4 ints
      real*8 vals(lenscr)
      integer labs(lenscr),bins(lenscr),ijv(lenbuf),klv(lenbuf)
      integer labels(2,lenbuf),iork(lenbuf),jorl(lenbuf)
c
      data mask10 /1023/
      save mask10
      parameter (debug=.false.)
c
c
      nnpsq=nnp**2
c
c     ----- initialise the sort routines -----
c
      if(debug) then
         write(iout,*) 'calling sorter:nnp,nnpsq ',nnp,nnpsq
      endif
      call sorter('start',acore,acore,maxcor,nnpsq,0,0,0,0,
     #           'sorted ao derivative integrals','scr',prnt)
c
c
      tot=0
      ptr=last
    1 continue
         call iosys('read real "unsorted derivative integrals" '//
     $              'from rdints',lenbuf,labels,ptr,' ')
         call iosys('read real "unsorted derivative integrals" '//
     $              'from rdints without rewinding',lenbuf,values,0,' ')
c
         n=abs(values(1))
         ptr=labels(1,1)
c
         do 2 i=2,n
            iork(i)=shiftr(labels(1,i),10)
            jorl(i)=and(labels(1,i),mask10)
    2    continue
         do 3 i=2,n
            ijv(i)=max(iork(i),jorl(i))
            jorl(i)=min(iork(i),jorl(i))
    3    continue
         do 4 i=2,n
            ijv(i)=ijv(i)*(ijv(i)-1)/2+jorl(i)
    4    continue
         do 5 i=2,n
            iork(i)=shiftr(labels(2,i),10)
            jorl(i)=and(labels(2,i),mask10)
    5    continue
         do 6 i=2,n
            klv(i)=max(iork(i),jorl(i))
            jorl(i)=min(iork(i),jorl(i))
    6    continue
         do 7 i=2,n
            klv(i)=klv(i)*(klv(i)-1)/2+jorl(i)
    7    continue
         do 23 i=1,n-1
            labs(i)=(ijv(i+1)-1)*nnp+klv(i+1)
            vals(i)=values(i+1)
   23    continue
         junk=n-1
         do 24 i=1,n-1
            labs(i+junk)=(klv(i+1)-1)*nnp+ijv(i+1)
            vals(i+junk)=values(i+1)
   24    continue
c
         count=2*n-2
         tot=tot+count
         if(tot.gt.3*nnpsq) then
            write(iout,*) 'warning: tot gt 3nnpsq ',tot,nnpsq
         endif
         if(tot.gt.6*nnpsq) then
            write(iout,*) ' problem tot gt 6nnpsq ',tot,nnpsq
            call lnkerr('m824:sort32')
         endif
         call sorter('with bin',acore,acore,0,count,labs,bins,vals,
     #          0,0,0,prnt)
c
c
      if (ptr.ge.0) go to 1
c
c     ----- finish the sort, ending on ints -----
c
      call sorter('end',acore,acore,0,0,0,0,0,0,0,0,prnt)
      if(debug) then
         write(iout,*) ' total words processed ',tot
      endif
c
c
      return
      end
