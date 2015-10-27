*deck @(#)chorre.f	5.1  11/6/94
      subroutine chorre(joutfg,korboc)
c
c  check the configurations to see if they satisfy the orbital
c  restrictions
c
      common /c1/ nbf,maxhsh,maxorb,ni,n,nbfp1,ibdcfg,icfgpn
      common /datarr/ kset(10),lmset(2,10),mmset(100)
      dimension joutfg(nbfp1,*)
c
      ibdcfg = 0
      lll = 0
      do 20 jj = 1,korboc
         norbs = kset(jj)
         neln = 0
         do 10 ki = 1,norbs
            lll = lll + 1
   10       neln = neln + joutfg(mmset(lll),n)
         if(neln.lt.lmset(1,jj) .or. neln.gt.lmset(2,jj))  go to 30
   20 continue
      return
   30 ibdcfg = 1
      return
      end
