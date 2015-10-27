*deck @(#)sorter1.f	5.1  11/6/94
      subroutine sorter1(n,labels,bins,values,lbin,vbin,binsiz,perbin,
     #                   offset,bintot,sortq)
c
c***begin prologue     sorter1
c***date written       850601   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           sorter dependent routines
c
c***author             saxe, paul,    (lanl)
c***source             @(#)sorter1.f	5.1   11/6/94
c***purpose            dependent routine of sorter used to place values
c                         in bins, writing out full bins.
c***description        #
c
c
c***references         (none)
c
c***routines called    iosys   (io)
c
c   common blocks:     (none)
c
c***end prologue       sorter1
c
      implicit integer (a-z)
c
      character*(*) sortq,sort*16
      real*8 values(n),vbin(*)
      integer labels(n),bins(n),lbin(*)
c
c     ----- bins are treated as linear arrays, so find offset -----
c
      junk=offset+1
c
cdir$ fastmd
c
      do 1 i=1,n
         bins(i)=(labels(i)-junk)/perbin*binsiz
    1 continue
c
cdir$ slowmd
c
c
c     ----- put each element in the proper bin -----
c
      do 2 i=1,n
         bin=bins(i)
         pt=lbin(1+bin)+1
c
c        ----- check if the bin is full -----
c
         if (pt.gt.binsiz) then
            sort=sortq
            call iosys('get write pointer of bins on '//sort,chain,0,0,
     #                                                             ' ')
            call iosys('write integer bins on '//sort//
     #                  ' without rewinding',binsiz,lbin(1+bin),0,' ')
            call iosys('write integer bins on '//sort//
     #                  ' without rewinding',wptoin(binsiz),
     #                  vbin(1+bin),0,' ')
            lbin(2+bin)=chain
            bintot=bintot+1
            pt=3
         end if
c
         lbin(1+bin)=pt
         lbin(pt+bin)=labels(i)
         vbin(pt+bin)=values(i)
    2 continue
c
c
      return
      end
