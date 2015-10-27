c
c     $Header: /msrc/proj/mss/tcgmsg/ipcv4.0/hpuxargs.f,v 1.1 1994/02/23 16:54:23 d3g681 Exp $
c
      integer function hpargc()
      hpargc = iargc() + 1
      end
      integer function hpargv(index, arg, maxlen)
      character*256 arg
      hpargv = igetarg(index,arg,maxlen)
      end
