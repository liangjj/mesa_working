*deck @(#)nxtask.f	5.1  11/28/95
      integer function nxtask(nproc)
      parameter (ichunk = 10)
      save icount, nleft
      data nleft, icount /0, 0/
c
c     wrapper round nxtval() to increase granularity
c     and thus reduce no. of requests to shared counter
c
      if(nproc.gt.0) then
         if(nleft.eq.0) then
            icount = nxtval(nproc) * ichunk
            nleft = ichunk
         endif
         nxtask = icount
         icount = icount + 1
         nleft = nleft -1
      else
          nleft = 0
          nxtask = 0
          junk = nxtval(nproc)
      endif
c
c     following does dumb static load balancing
c
c$$$      if(nproc.gt.0) then
c$$$         if (nleft .eq. 0) then
c$$$            icount = nodeid()
c$$$            nleft = 1
c$$$         endif
c$$$         nxtask = icount
c$$$         icount = icount + nnodes()
c$$$      else
c$$$          nleft = 0
c$$$          nxtask = 0
c$$$      endif
      end
