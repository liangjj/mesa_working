*deck %W%  %G%
      subroutine rd39(vec,nbf,nob,nset,itap39)
C
C***Begin prologue
C***Date written       871022   (yymmdd)
C***Revision date      yymmdd   (yymmdd)
C
C***Keywords
C***Author             Lengsfield, Byron (BRL)
C***Source             %W%   %G%
C
C***Purpose
C
C***Description
C
C***References
C
C***Routines called    (none)
C
C***End prologue
C
      implicit real*8 (a-h,o-z)
      dimension vec(nbf,2)
c
c   ipt39 is common to programs wr39 and rd39
c
cc      common / ipt39 / ipoint(4,101)
c
      ntot=nbf*nob
c
      call lnkerr(' RD39 called ')
c
c.io      call sread(itap39,vec,intowp(ntot))
c
      return
      end
