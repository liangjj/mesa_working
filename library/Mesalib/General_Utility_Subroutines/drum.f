*deck @(#)drum.f	5.1  11/6/94
      subroutine drum
c***begin prologue     drum
c***date written       850601  (yymmdd)
c***revision date      870207  (yymmdd)
c
c    7 february 1987   pws at lanl
c       extending 'names' to 128 characters to handle unix file names
c
c    1 december 1986   pws at lanl
c       changing 'names' to a character array.
c
c***keywords           link, external files, internal files, input/output
c***author             martin, richard (lanl)
c***source
c***purpose            starts up a new link.
c***description
c                      call drum
c                      the sequence of events:
c                        lxopen  opens the external files, (inp,out).
c                        liopen  opens and initializes the rwf file,
c                                restores the global common area,
c                                processes file name replacement and memory
c                                requests.
c
c***references
c***routines called    lxopen(mdutil), liopen(util)
c***end prologue       drum
c
      implicit integer(a-z)
      common/io/inp,iout
      character*128 names(50)
c
c
      call lxopen(names,mchsiz,maxsiz)
      call liopen(names,maxsiz)
c
c     we may want to process system dump options, assembler'ok',etc.
c     here someday.
c
c
      return
      end
