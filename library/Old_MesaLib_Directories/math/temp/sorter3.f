*deck @(#)sorter3.f	5.1  11/6/94
      subroutine sorter3(nbins,binsiz,chain,lbin,vbin,sorted,perbin,
     #                   accume,file,unit,ndata,sortq,global)
c
c***begin prologue     sorter3
c***date written       850601   (yymmdd)
c***revision date      910917   (yymmdd)
c   17 september 1991  rlm at lanl
c      "inlining" the call to scatter. some machines have trouble with
c      the negative addresses which can sometimes be generated via 'offset'.
c***keywords           sorter dependent routines
c
c***author             saxe, paul,    (lanl)
c***source             @(#)sorter3.f	5.1   11/6/94
c***purpose            sorter dependent routine to chain back through
c                         sort bins, constructing the final file.
c
c***description        #
c
c
c***references         (none)
c
c***routines called    lnkerr  (mdutil)
c                      iosys   (io)
c
c   common blocks:     io
c
c***end prologue       sorter3
c
      implicit integer (a-z)
c
      real*8 vbin(binsiz),sorted(perbin)
      integer chain(nbins),lbin(binsiz)
      character*32 file
      character*16 unit,sort,sortq*(*)
c
      common /io/     inp,iout
c
c     ----- rewind or create output file -----
c
      if (accume.ge.0) then
          call iosys('get maximum length of "'//file//'" on '//unit,
     #                                          len,0,0,' ')
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
      call iosys('rewind "'//file//'" on '//unit,0,0,0,' ')
c
c     ----- loop through each type of bin -----
c
      sort=sortq
      do 10 bin=1,nbins
c
         offset=(bin-1)*perbin+global
c
         if (accume.ge.0) then
            call rzero(sorted,perbin)
         else
            if (bin.lt.nbins) then
               junk=perbin
            else
               junk=ndata-(nbins-1)*perbin
            end if
            call iosys('read real "'//file//'" from '//unit,
     #                 junk,sorted,offset,' ')
         end if
c
c        ----- starting with last bin, chain back through bins -----
c
         locbin=chain(bin)
         if (locbin.lt.0) go to 5
    1    continue
            call iosys('read integer bins from '//sort,binsiz,lbin,
     #                 locbin,' ')
            call iosys('read integer bins from '//sort//
     #                 ' without rewinding',wptoin(binsiz),vbin,0,' ')
            n=lbin(1)
            locbin=lbin(2)
            if (accume.ge.0) then
               do 2 i=3,n
                  sorted(lbin(i)-offset)=vbin(i)
    2          continue
            else
               do 3 i=3,n
                  sorted(lbin(i)-offset)=sorted(lbin(i)-offset)+
     #                                     vbin(i)
    3          continue
            end if
         if (locbin.ge.0) go to 1
c
    5    continue
c
c        ----- write this part of the output out -----
c
         if (bin.lt.nbins) then
            junk=perbin
         else
            junk=ndata-(nbins-1)*perbin
         end if
c
c         write(iout,*) file, junk, offset
c         call lnkerr('quit')
         call iosys('write real "'//file//'" on '//unit,junk,
     #               sorted,offset,' ')
c
   10 continue
c
c
      return
      end
