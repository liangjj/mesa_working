*deck @(#)sort64.f	5.1  11/6/94
      subroutine sort64(labels,values,acore,lenbuf,nnp,ints,iout,
     $                 maxcor,vals,labs,bins,lenscr,ijv,klv,
     $                 iork,jorl,prnt,ops,killr)
c***begin prologue     sort64.f
c***date written       yymmdd  
c***revision date      11/6/94      
c
c***keywords           
c***author             saxe, paul(lanl)
c***source             @(#)sort64.f	5.1   11/6/94
c***purpose            64-bit integral sort 
c***description
c   19 november 1992    rlm at lanl
c      removing the input data file after the data has been sorted onto
c      the scratch file. this is an attempt to cut down on the disk
c      requirements of the sort.
c   15 february 1987    pws at lanl
c      adding an option to print the integrals with labels.
c
c   nb.  iork and jorl are equivalenced to labs by the calling sequence
c
c***references
c
c***routines called
c
c***end prologue       sort64.f
      implicit none
c     --- input variables -----
      integer lenscr,lenbuf,nnp,iout,maxcor
      character*(*) ops
      character*4 ints
      logical prnt, killr 
c     --- input arrays (unmodified) ---
c     --- input arrays (scratch) ---
      integer labs(lenscr),bins(lenscr),ijv(lenbuf),klv(lenbuf)
      integer iork(lenbuf),jorl(lenbuf)
      integer labels(lenbuf)
      real*8 values(lenbuf),acore(maxcor)
      real*8 vals(lenscr)
c     --- output arrays ---
c     --- output variables ---
c     --- scratch arrays ---
c     --- local variables ---
      integer mask10
      integer and
      integer wpadti,n,i,shiftr,ilbl,jlbl,klbl,llbl,count,junk
      logical logkey
c
c     data mask10 /o'1777'/
      parameter (mask10=1023)
c
c     --- initialise the sort routines ---
      call sorter('start',acore,acore,wpadti(maxcor),nnp**2,0,0,0,0,
     $            'sorted ao integrals',ints,prnt)
c
c
    1 continue
         call iosys('read real "unsorted ao integrals" from rints '//
     $              'without rewinding',lenbuf,labels,0,' ')
         call iosys('read real "unsorted ao integrals" from rints '//
     $              'without rewinding',lenbuf,values,0,' ')
c
         n=abs(values(1))
c
         do 2 i=2,n
            iork(i)=shiftr(labels(i),30)
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
         if (logkey(ops,'print=m330=ints',.false.,' ')) then
            do 801 i=1,n-1
               ilbl=and(shiftr(labels(i+1),30),mask10)
               jlbl=and(shiftr(labels(i+1),20),mask10)
               klbl=and(shiftr(labels(i+1),10),mask10)
               llbl=and(labels(i+1),mask10)
               write (iout,800) i,ilbl,jlbl,klbl,llbl,values(i+1)
  800          format (1x,i5,3x,4i3,f20.9)
  801       continue
         end if
c
         count=2*n-2
         call sorter('with bin',acore,acore,0,count,labs,bins,vals,
     $          0,0,0,prnt)
c
c
      if (values(1).gt.0.0d+00) go to 1
c
c        --- finish the sort, ending on ints ---
c             destroy the input file to free disk space for the output file.
      if(killr) then
         call iosys('destroy rints',0,0,0,' ')
      endif
      call sorter('end',acore,acore,0,0,0,0,0,0,0,0,prnt)
c
c
      return
      end
