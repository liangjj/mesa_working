*deck @(#)pm2001.f	5.1  11/6/94
      program m2001
c***begin prologue     pm2001.f
c***date written       850601   (yymmdd)
c***revision date      11/6/94
c***keywords           m2001, link 2001, archive, exit
c***author             martin, richard (lanl)
c***source             @(#)pm2001.f	5.1   11/6/94
c***purpose            the termination link for mesa.  it makes sure the
c                      chk file is in order, and may prepare an archive
c                      entry.
c***description
c     this is the last link executed in a mesa run.
c     it's primary purpose is to make sure the appropriate information
c     about the run is saved in the chk file to permit restarts,
c     subsequent interrogations, etc.
c     it may also do all, or none, of the following:
c        prepare and print an archive entry.
c        store the archive entry.
c        print a timing summary for all links previously executed.
c
c***references
c
c***routines called
c***end prologue       m2001
      implicit none
c     --- input variables -----
c     --- input arrays (unmodified) ---
c     --- input arrays (scratch) ---
c     --- output arrays ---
c     --- output variables ---
c     --- scratch arrays ---
c     --- local variables ---
      integer inp,iout
      character*4096 ops
c
      common/io/inp,iout
      call drum
c
c     --- retrieve the link options.
      call iosys('read character options from rwf',-1,0,0,ops)
c
c     --- checkpoint the run.
      call chkpnt('save')
c
c     --- prepare the archive entry.
      call archiv
c
c     --- print the global time summary.
      call tsumry
c
c     --- exit.
      call chainx(0)
c
c
      end
