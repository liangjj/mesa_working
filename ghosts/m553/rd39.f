      subroutine rd39(vec,nbf,nob,nset,itap39)
c
c***begin prologue
c***date written       871022   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c
c***keywords
c***author             lengsfield, byron (brl)
c***source             %W%   %G%
c
c***purpose
c
c***description
c
c***references
c
c***routines called    (none)
c
c***end prologue
c
      dimension vec(nbf,2)
c
c   ipt39 is common to programs wr39 and rd39
c
cc      common / ipt39 / ipoint(4,101)
c
      ntot=nbf*nob
c
      call lnkerr(' rd39 called ')
c
c.io      call sread(itap39,vec,intowp(ntot))
c
      return
      end
