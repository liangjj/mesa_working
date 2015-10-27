      subroutine wr39(vec,nbf,nob,nset,itap39)
c
      dimension vec(nbf,2)
c
c   ipt39 is common to programs wr39 and rd39
c
cc      common / ipt39 / ipoint(4,101)
c
      ntot=nbf*nob
      call lnkerr(' wr39 called ')
c
c
      return
      end
