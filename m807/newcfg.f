*deck @(#)newcfg.f	5.1  11/6/94
      function newcfg(j,nhash)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c  convert a list of 0,1,2's to a unique ternary number; calculate a
c  certain modulus as an array index for checking the array for the
c  ternary number; if not encountered, enter it in the hash table
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
*mdc*if harris
*      integer*6 nhash, num1u, num2u, num3u
*      parameter (i1u=28, i2l=29, i2u=56, i3l=57, i3u=84)
*mdc*elseif cray
*      parameter (i1u=38, i2l=39, i2u=76, i3l=77, i3u=114)
*mdc*else
      parameter (i1u=18, i2l=19, i2u=36, i3l=37, i3u=54)
*mdc*endif
c
c mesa
c
      common /io/inp,iout
c
      common /c1/ nbf,maxhsh,maxorb,ni,n,nbfp1,ibdcfg,icfgpn
      dimension j(*), nhash(*)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     first convert j to a ternary number
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      kflag = 0
      if( nbf .gt. i1u )   go to 140
      num1 = 0
      do 110 ibf = 1,nbf
  110 num1 = num1*3 + j(ibf)
      ndx = mod(num1,24576) + 1
  120 if( nhash(ndx) .ne. 0 )   go to 130
      newcfg = 0
      nhash(ndx) = num1
      return
  130 if( nhash(ndx) .eq. num1 )  go to 270
      ndx = ndx + 1
      if( ndx .le. maxhsh )   go to 120
      if( kflag .ne. 0 )   go to 280
      ndx = 1
      kflag = 1
      go to 120
  140 if( nbf .gt. i2u )   go to 210
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     calculate num1, num2 from max of i1u entries each
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      num1 = 0
      do 150 ibf = 1,i1u
  150 num1 = num1*3 + j(ibf)
      num2 = 0
      do 160 ibf = i2l,nbf
  160 num2 = num2*3 + j(ibf)
      ndx = 2*mod(num1+num2,12288) + 1
  170 if( nhash(ndx) .eq. 0 )   go to 180
      if(nhash(ndx).eq.num1 .and. nhash(ndx+1).eq.num2)  go to 270
      ndx = ndx + 2
      if( ndx .lt. maxhsh )   go to 170
      if( kflag .ne. 0 )   go to 280
      ndx = 1
      kflag = 1
      go to 170
  180 nhash(ndx)   = num1
      nhash(ndx+1) = num2
      newcfg = 0
      return
  210 if( nbf .gt. i3u )   go to 290
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     calculate num1, num2, num3 from max of i1u entries each
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      num1 = 0
      do 220 ibf = 1,i1u
  220 num1 = num1*3 + j(ibf)
      num2 = 0
      do 230 ibf = i2l,i2u
  230 num2 = num2*3 + j(ibf)
      num3 = 0
      do 240 ibf = i3l,nbf
  240 num3 = num3*3 + j(ibf)
      ndx = 3*mod(num1+num2+num3,8192) + 1
  250 if( nhash(ndx) .eq. 0 )   go to 260
      if(nhash(ndx)  .eq.num1 .and. nhash(ndx+1).eq.num2 .and.
     1   nhash(ndx+2).eq.num3)  go to 270
      ndx = ndx + 3
      if( ndx .lt. maxhsh )   go to 250
      if( kflag .ne. 0 )   go to 280
      ndx = 1
      kflag = 1
      go to 250
  260 nhash(ndx)   = num1
      nhash(ndx+1) = num2
      nhash(ndx+2) = num3
      newcfg = 0
      return
  270 newcfg = 1
      return
  280 write (iout,310)
      stop
  290 write (iout,300)   nbf, i3u
      stop
  300 format(/11x,'number of basis functions for newcfg restricted',
     $   ' to',i4,' nbf =',i7)
  310 format(///11x,'hash table for configurations is full.')
      end
