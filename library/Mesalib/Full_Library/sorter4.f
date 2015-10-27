*deck @(#)sorter4.f	5.1  11/6/94
      subroutine sorter4(vals,vbin,val2,sorted,labs,lbin,bins,lab2,
     #                   chain,chain2,binsz1,binsz2,nbins,perbn1,
     #                   perbn2,accume,file,unit,ndata)
c
c***begin prologue     sorter4
c***date written       850601   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           sorter dependent routines
c
c***author             saxe, paul,    (lanl)
c***source             @(#)sorter4.f	5.1   11/6/94
c***purpose            sorter routine to complete a two-pass sort by
c                         sorting each of the first pass bins and then
c                         formaing the final file.
c***description        #
c
c
c***references         (none)
c
c***routines called    lnkerr  (mdutil)
c                      iosys   (io)
c                      sorter1 (math)
c                      sorter2 (math)
c                      sorter3 (math)
c
c   common blocks:     (none)
c
c***end prologue       sorter4
c
      implicit integer (a-z)
c
      character*32 file
      character*16 unit
      real*8 vals(binsz1),vbin(binsz2,nbins),val2(binsz2)
      integer labs(binsz1),lbin(binsz2,nbins),chain(nbins)
      integer lab2(binsz2),chain2(nbins),bins(binsz1)
c
c     ----- rewind or create output file -----
c
      if (accume.ge.0) then
          call iosys('get maximum length of "'//file//'" on '//unit,
     #                                 len,0,0,' ')
         if (len.le.0) then
            call iosys('create real "'//file//'" on '//unit,
     $           ndata,0,0,' ')
         else
            if (len.lt.ndata) then
               call lnkerr(file//' on '//unit//' already exists '//
     #                     'but is not long enough')
            end if
         end if
      end if
c
c
c     ----- loop through the bin types from first phase -----
c
      call iosys('rewind bins on sort',0,0,0,' ')
c
      do 10 bin=1,nbins
c
c        ----- initialize the phase two bins -----
c
         do 1 i=1,nbins
            lbin(1,i)=2
            lbin(2,i)=-1
    1    continue
c
         offset=(bin-1)*perbn1
c
c        ----- starting with the last bin, chain back through bins ----
c
         locbin=chain(bin)
c
         if (locbin.lt.0) go to 5
c
      call iosys('rewind bins on sort1',0,0,0,' ')
c
    2    continue
            call iosys('read integer bins from sort',binsz1,labs,
     #                 locbin,' ')
            call iosys('read integer bins from sort without rewinding',
     #                  wptoin(binsz1),vals,0,' ')
            n=labs(1)
            locbin=labs(2)
            call sorter1(n-2,labs(3),bins,vals(3),lbin,vbin,binsz2,
     #                   perbn2,offset,bintot,'sort1')
         if (locbin.ge.0) go to 2
c
    5    continue
c
c        ----- flush the final bins -----
c
         call sorter2(lbin,vbin,binsz2,nbins,chain2,bintot,'sort1')
c
         junk=min(perbn1,ndata-offset)
         call sorter3(nbins,binsz2,chain2,lab2,val2,sorted,perbn2,
     #                accume,file,unit,junk,'sort1',offset)
c
   10 continue
c
c
      return
      end
