*deck @(#)getkep.f	1.1  7/9/91
      subroutine getkep(line,keep,nkeep)
      implicit integer(a-z)
c
      character*(*) line
      character*3 keep(10)
c
      logical logkey
c
      nkeep=0
      if(logkey(line,'keep=s',.false.,' ')) then
         nkeep=nkeep+1
         keep(nkeep)='s'
      end if
      if(logkey(line,'keep=x',.false.,' ')) then
         nkeep=nkeep+1
         keep(nkeep)='x'
      end if
      if(logkey(line,'keep=y',.false.,' ')) then
         nkeep=nkeep+1
         keep(nkeep)='y'
      end if
      if(logkey(line,'keep=z',.false.,' ')) then
         nkeep=nkeep+1
         keep(nkeep)='z'
      end if
      if(logkey(line,'keep=xx',.false.,' ')) then
         nkeep=nkeep+1
         keep(nkeep)='xx'
      end if
      if(logkey(line,'keep=yy',.false.,' ')) then
         nkeep=nkeep+1
         keep(nkeep)='yy'
      end if
      if(logkey(line,'keep=zz',.false.,' ')) then
         nkeep=nkeep+1
         keep(nkeep)='zz'
      end if
      if(logkey(line,'keep=xy',.false.,' ')) then
         nkeep=nkeep+1
         keep(nkeep)='xy'
      end if
      if(logkey(line,'keep=xz',.false.,' ')) then
         nkeep=nkeep+1
         keep(nkeep)='xz'
      end if
      if(logkey(line,'keep=yz',.false.,' ')) then
         nkeep=nkeep+1
         keep(nkeep)='yz'
      end if
c..bhl sept 27
      if(logkey(line,'keep=xxx',.false.,' ')) then
         nkeep=nkeep+1
         keep(nkeep)='xxx'
      end if
      if(logkey(line,'keep=yyy',.false.,' ')) then
         nkeep=nkeep+1
         keep(nkeep)='yyy'
      end if
      if(logkey(line,'keep=zzz',.false.,' ')) then
         nkeep=nkeep+1
         keep(nkeep)='zzz'
      end if
      if(logkey(line,'keep=xxz',.false.,' ')) then
         nkeep=nkeep+1
         keep(nkeep)='xxz'
      end if
      if(logkey(line,'keep=xyy',.false.,' ')) then
         nkeep=nkeep+1
         keep(nkeep)='xyy'
      end if
      if(logkey(line,'keep=xyz',.false.,' ')) then
         nkeep=nkeep+1
         keep(nkeep)='xyz'
      end if
      if(logkey(line,'keep=xzz',.false.,' ')) then
         nkeep=nkeep+1
         keep(nkeep)='xzz'
      end if
      if(logkey(line,'keep=yyz',.false.,' ')) then
         nkeep=nkeep+1
         keep(nkeep)='yyz'
      end if
      if(logkey(line,'keep=yzz',.false.,' ')) then
         nkeep=nkeep+1
         keep(nkeep)='yzz'
      end if
      if(logkey(line,'keep=xxy',.false.,' ')) then
         nkeep=nkeep+1
         keep(nkeep)='xxy'
      end if
c
      return
      end
