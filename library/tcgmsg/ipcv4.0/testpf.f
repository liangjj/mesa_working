c     $Header: /msrc/proj/mss/tcgmsg/ipcv4.0/testpf.f,v 1.1 1994/02/23 16:55:29 d3g681 Exp $
      character*60 fname
      call pbeginf
      fname = ' '
      write(fname,'(a,i3.3)') '/tmp/pfcopy.test',nodeid()
      call pfcopy(5, 0, fname)
      call pend
      end
