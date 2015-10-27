*deck @(#)nlmxyz.f	1.1 9/7/91
c***begin prologue     nlmxyz
c***date written       890404   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           nlmxyz, link 6001, nlm 
c***author             schneider, barry (lanl)
c***source             m6001
c***purpose            nlm to xyz of cartesian gaussians through f orbitals
c***description        fills up array with character values of cartesian
c***                   aos up to f functions bases on their n, l, m
c*** 
c
c***references       
c
c***routines called    iosys, util and mdutil
c***end prologue       nlmxyz
      subroutine nlmxyz
      common /ctype/ ctype(0:3,0:3,0:3)
      character *3 ctype
      ctype(0,0,0)='s'
      ctype(1,0,0)='x'
      ctype(0,1,0)='y'
      ctype(0,0,1)='z'
      ctype(2,0,0)='xx'
      ctype(0,2,0)='yy'
      ctype(0,0,2)='zz'
      ctype(1,1,0)='xy'
      ctype(1,0,1)='xz'
      ctype(0,1,1)='yz'
      ctype(3,0,0)='xxx'
      ctype(0,3,0)='yyy'
      ctype(0,0,3)='zzz'
      ctype(2,1,0)='xxy'
      ctype(2,0,1)='xxz'
      ctype(1,2,0)='xyy'
      ctype(0,2,1)='yyz'
      ctype(1,0,2)='xzz'
      ctype(0,1,2)='yzz'
      ctype(1,1,1)='xyz'
      return
      end
