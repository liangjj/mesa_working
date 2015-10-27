*deck @(#)getdrp.f	1.1  7/9/91
      subroutine getdrp(line,drop,nzero)
      implicit integer(a-z)
c
      character*(*) line
      character*3 drop(10)
c
      logical logkey
c
      nzero=0
      if(logkey(line,'drop=s',.false.,' ')) then
         nzero=nzero+1
         drop(nzero)='s'
      end if
      if(logkey(line,'drop=x',.false.,' ')) then
         nzero=nzero+1
         drop(nzero)='x'
      end if
      if(logkey(line,'drop=y',.false.,' ')) then
         nzero=nzero+1
         drop(nzero)='y'
      end if
      if(logkey(line,'drop=z',.false.,' ')) then
         nzero=nzero+1
         drop(nzero)='z'
      end if
      if(logkey(line,'drop=xx',.false.,' ')) then
         nzero=nzero+1
         drop(nzero)='xx'
      end if
      if(logkey(line,'drop=yy',.false.,' ')) then
         nzero=nzero+1
         drop(nzero)='yy'
      end if
      if(logkey(line,'drop=zz',.false.,' ')) then
         nzero=nzero+1
         drop(nzero)='zz'
      end if
      if(logkey(line,'drop=xy',.false.,' ')) then
         nzero=nzero+1
         drop(nzero)='xy'
      end if
      if(logkey(line,'drop=xz',.false.,' ')) then
         nzero=nzero+1
         drop(nzero)='xz'
      end if
      if(logkey(line,'drop=yz',.false.,' ')) then
         nzero=nzero+1
         drop(nzero)='yz'
      end if
c..bhl sept 27
      if(logkey(line,'drop=xxx',.false.,' ')) then
         nzero=nzero+1
         drop(nzero)='xxx'
      end if
      if(logkey(line,'drop=yyy',.false.,' ')) then
         nzero=nzero+1
         drop(nzero)='yyy'
      end if
      if(logkey(line,'drop=zzz',.false.,' ')) then
         nzero=nzero+1
         drop(nzero)='zzz'
      end if
      if(logkey(line,'drop=xxz',.false.,' ')) then
         nzero=nzero+1
         drop(nzero)='xxz'
      end if
      if(logkey(line,'drop=xyy',.false.,' ')) then
         nzero=nzero+1
         drop(nzero)='xyy'
      end if
      if(logkey(line,'drop=xyz',.false.,' ')) then
         nzero=nzero+1
         drop(nzero)='xyz'
      end if
      if(logkey(line,'drop=xzz',.false.,' ')) then
         nzero=nzero+1
         drop(nzero)='xzz'
      end if
      if(logkey(line,'drop=yyz',.false.,' ')) then
         nzero=nzero+1
         drop(nzero)='yyz'
      end if
      if(logkey(line,'drop=yzz',.false.,' ')) then
         nzero=nzero+1
         drop(nzero)='yzz'
      end if
      if(logkey(line,'drop=xxy',.false.,' ')) then
         nzero=nzero+1
         drop(nzero)='xxy'
      end if
c
      return
      end
