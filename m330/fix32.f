*deck  @(#)fix32.f	5.1 11/6/94
      subroutine fix32(labels,values,acore,lenbuf,nnp,ints,iout,
     $                 maxcor,vals,labs,bins,lenscr,ijv,klv,
     $                 iork,jorl,prnt,ops,pkindx,killr)
c***begin prologue     fix32.f
c***date written       yymmdd  
c***revision date      11/6/94      
c
c***keywords           
c***author             saxe, paul(lanl)
c***source             @(#)fix32.f	5.1   11/6/94
c***purpose            32-bit integral sort
c***description
c
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
c***end prologue       fix32.f
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
      integer labels(2,lenbuf),pkindx(*)
      real*8 values(lenbuf),acore(maxcor)
      real*8 vals(lenscr)
c     --- output arrays ---
c     --- output variables ---
c     --- scratch arrays ---
c     --- local variables ---
      integer mask10
      integer and
      integer wpadti,n,i,shiftr,icnt,ilbl,jlbl,klbl,llbl,count,junk
      logical logkey
c
c     data mask10 /o'1777'/
      parameter (mask10=1023)
c
c     --- initialise the sort routines ---
      call sorter('start',acore,acore,wpadti(maxcor),nnp**2,0,0,0,0,
     $           'sorted ao integrals',ints,prnt)
c
    1 continue
         call iosys('read real "unsorted ao integrals" from rints '//
     $              'without rewinding',lenbuf,labels,0,' ')
         call iosys('read real "unsorted ao integrals" from rints '//
     $              'without rewinding',lenbuf,values,0,' ')
c
         n=abs(values(1))
c
c        --- pack the incoming integrals ---
c        --- pull out labels ---
         do 11 i=2,n
            iork(i)=shiftr(labels(1,i),10)
            jorl(i)=and(labels(1,i),mask10)
            ijv(i)=shiftr(labels(2,i),10)
            klv(i)=and(labels(2,i),mask10)
   11    continue
c
c        --- check if integral is to be dropped ---
         do 12 i=2,n
            iork(i)=pkindx(iork(i))*pkindx(jorl(i))*
     $              pkindx(ijv(i))*pkindx(klv(i))
   12    continue
c
c        --- if not put old label in iork and values --- 
         icnt=2
         do 13 i=2,n
            if(iork(i).ne.0) then
               labels(1,icnt)=labels(1,i)
               labels(2,icnt)=labels(2,i)
               values(icnt)=values(i)
               icnt=icnt+1
            end if
   13    continue
c
         n=icnt-1
         if(n.ne.1) then
c
c           --- retrieve old indices ---
            do 2 i=2,n
               iork(i)=shiftr(labels(1,i),10)
               jorl(i)=and(labels(1,i),mask10)
    2       continue
c
c           --- make new indices ---
            do 21 i=2,n
               iork(i)=pkindx(iork(i))
               jorl(i)=pkindx(jorl(i))
   21       continue
c
            do 3 i=2,n
               ijv(i)=max(iork(i),jorl(i))
               jorl(i)=min(iork(i),jorl(i))
    3       continue
            do 4 i=2,n
               ijv(i)=ijv(i)*(ijv(i)-1)/2+jorl(i)
    4       continue
c
c           --- now for the other two indices ---
            do 5 i=2,n
               iork(i)=shiftr(labels(2,i),10)
               jorl(i)=and(labels(2,i),mask10)
    5       continue
c
            do 51 i=2,n
               iork(i)=pkindx(iork(i))
               jorl(i)=pkindx(jorl(i))
   51       continue
c
            do 6 i=2,n
              klv(i)=max(iork(i),jorl(i))
              jorl(i)=min(iork(i),jorl(i))
    6       continue
            do 7 i=2,n
               klv(i)=klv(i)*(klv(i)-1)/2+jorl(i)
    7       continue
c
c           --- put in scratch buffer ---
            do 23 i=1,n-1
               labs(i)=(ijv(i+1)-1)*nnp+klv(i+1)
               vals(i)=values(i+1)
   23       continue
            junk=n-1
            do 24 i=1,n-1
               labs(i+junk)=(klv(i+1)-1)*nnp+ijv(i+1)
               vals(i+junk)=values(i+1)
   24       continue
c
            if (logkey(ops,'print=m330=ints',.false.,' ')) then
               do 801 i=1,n-1
                  ilbl=shiftr(labels(1,i+1),10)
                  jlbl=and(labels(1,i+1),mask10)
                  klbl=shiftr(labels(2,i+1),10)
                  llbl=and(labels(2,i+1),mask10)
                  write (iout,800) i,ilbl,jlbl,klbl,llbl,values(i+1)
  800             format (1x,i5,3x,4i3,f20.9)
  801          continue
            end if
c
c           --- away to sorter ---     
            count=2*n-2
            call sorter('with bin',acore,acore,0,count,labs,bins,vals,
     $                   0,0,0,prnt)
c
         endif
c
      if (values(1).gt.0.0d+00) go to 1
c
c     --- finish the sort, ending on ints ---
c         destroy input file to free disk space for output file.
      if(killr) then
         call iosys('destroy rints',0,0,0,' ')
      endif
      call sorter('end',acore,acore,0,0,0,0,0,0,0,0,prnt)
c
c
      return
      end