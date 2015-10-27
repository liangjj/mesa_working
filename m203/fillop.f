*deck @(#)fillop.f	1.1  11/30/90
      subroutine fillop(ops,newop)
c***begin prologue     fillop
c***date written       850601  yymmdd
c***revision date      yymmdd  yymmdd
c***keywords           options
c***author             martin, richard (lanl)
c***source             m1
c***purpose            adds an option to the option string.
c***description
c     call fillop(ops,newop)
c       ops     the options string.
c       newop   the option to add.
c
c     this routine will add newop to the ops string if it doesn't
c     already exist there.
c***references
c***routines called    cskipb(chr), pakstr(chr)
c***end prologue       fillop
      implicit integer(a-z)
      character*(*) ops,newop
      character*1 comma,blank
c
      data comma/','/, blank/' '/
c
c     first find the end of the options string.
      lenops=cskipb(ops,blank)
c
c     pack down the new option.
      call pakstr(newop,lennew)
      if(lennew.eq.0) return
c
c     see if the new option already exists in the string.
      if(index(ops,newop(:lennew)).ne.0) return
c
c     place it in the options string followed by a comma.
      ops(lenops+1:)=newop(:lennew)//comma
c
c
      return
      end
