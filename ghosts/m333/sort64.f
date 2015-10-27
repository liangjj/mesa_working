*deck @(#)sort64.f	1.1  11/30/90
      subroutine sort64(labels,values,acore,lenbuf,num,nnp,ints,iout,
     #                 maxcor,vals,labs,bins,lenscr,ijv,klv,
     #                 iork,jorl,prnt,nder)
c
c   nb.  iork and jorl are equivalenced to labs by the calling sequence
c
      implicit integer (a-z)
c
      logical prnt
      real*8 values(lenbuf),acore(maxcor)
      real*8 temp
      character*4 ints
      real*8 vals(lenscr)
      integer labs(lenscr),bins(lenscr),ijv(lenbuf),klv(lenbuf)
      integer labels(lenbuf),iork(lenbuf),jorl(lenbuf)
c
      parameter (mask10=1023)
c     data mask10/1777b/
c
c     ----- start timing -----
c
      nnpsq=nnp**2
c
c     ----- initialise the sort routines -----
c
      call sorter('start',acore,acore,maxcor,nnpsq*nder,0,0,0,0,
     #           'sorted ao derivative integrals',ints,prnt)
c
c
    1 continue
         call iosys('read real "unsorted derivative integrals" '//
     $     'from ints without rewinding',lenbuf,labels,0,' ')
         call iosys('read real "unsorted derivative integrals" '//
     $        'from ints without rewinding',lenbuf,values,0,' ')
c
         n=abs(values(1))
c
c        shift1=shiftl(1,30)
c        shift2=shiftl(1,20)
         do 2 i=2,n
            iork(i)=and(shiftr(labels(i),30),mask10)
            jorl(i)=and(shiftr(labels(i),20),mask10)
    2    continue
         do 3 i=2,n
            ijv(i)=max(iork(i),jorl(i))
            jorl(i)=min(iork(i),jorl(i))
    3    continue
         do 4 i=2,n
            ijv(i)=ijv(i)*(ijv(i)-1)/2+jorl(i)
    4    continue
         do 5 i=2,n
            iork(i)=and(shiftr(labels(i),10),mask10)
            jorl(i)=and(labels(i),mask10)
    5    continue
         do 6 i=2,n
            klv(i)=max(iork(i),jorl(i))
            jorl(i)=min(iork(i),jorl(i))
    6    continue
         do 7 i=2,n
            klv(i)=klv(i)*(klv(i)-1)/2+jorl(i)
    7    continue
c        do 21 i=2,n
c           ijv(i)=labels(i)/shift1
c           ijv(i)=ijv(i)*(ijv(i)-1)/2+and(labels(i)/shift2,mask10)
c  21    continue
c        shift1=shiftl(1,10)
c        do 22 i=2,n
c           klv(i)=and(labels(i)/shift1,mask10)
c           klv(i)=klv(i)*(klv(i)-1)/2+and(labels(i),mask10)
c  22    continue
         do 9 i=2,n
            labels(i)=shiftr(labels(i),40)*nnpsq-nnpsq
    9    continue
         do 23 i=1,n-1
            labs(i)=(ijv(i+1)-1)*nnp+klv(i+1)+labels(i+1)
            vals(i)=values(i+1)
   23    continue
         junk=n-1
         do 24 i=1,n-1
            labs(i+junk)=(klv(i+1)-1)*nnp+ijv(i+1)+labels(i+1)
            vals(i+junk)=values(i+1)
   24    continue
c
         count=2*n-2
c        call sorter('with bin',acore,acore,0,0,0,0, 0,0,0,0,
c    #                vals,labs,bins,count,prnt)
         call sorter('with bin',acore,acore,0,count,labs,bins,vals,
     #          0,0,0,prnt)
c
c
      if (values(1).gt.0.0d+00) go to 1
c
c     ----- finish the sort, ending on ints -----
c
c     call sorter('end',acore,acore,0,0,0,0, 0,0,0,0, 0,0,0,0,prnt)
      call sorter('end',acore,acore,0,0,0,0,0,0,0,0,prnt)
c
c     ----- stop timing -----
c
c
c
      return
      end
