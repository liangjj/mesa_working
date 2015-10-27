*deck @(#)sorter2.f	5.1  11/6/94
      subroutine sorter2(lbin,vbin,binsiz,nbins,chain,bintot,sortq)
c
c***begin prologue     sorter2
c***date written       850601   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           sorter dependent routines
c
c***author             saxe, paul,    (lanl)
c***source             @(#)sorter2.f	5.1   11/6/94
c***purpose            sorter dependent routine to write the partially
c                         filled bins to disc at the end of an internal
c                         sort.
c***description        #
c
c
c***references         (none)
c
c***routines called    iosys   (io)
c
c   common blocks:     (none)
c
c***end prologue       sorter2
c
      implicit integer (a-z)
c
      character*(*) sortq,sort*16
      real*8 vbin(binsiz,nbins)
      integer lbin(binsiz,nbins),chain(nbins)
c
c     ----- write out the final partially filled bins, and put the
c           location of the last bin of each type in chain
c
      sort=sortq
      do 1 bin=1,nbins
         if (lbin(1,bin).le.2) then
            chain(bin)=lbin(2,bin)
         else
            call iosys('get write pointer of bins on '//sort,
     #                  chain(bin),0,0,' ')
            call iosys('write integer bins on '//sort//
     #                 ' without rewinding',binsiz,lbin(1,bin),0,' ')
            call iosys('write integer bins on '//sort//
     #                 ' without rewinding',wptoin(binsiz),
     #                 vbin(1,bin),0,' ')
            bintot=bintot+1
         end if
    1 continue
c
c
      return
      end
