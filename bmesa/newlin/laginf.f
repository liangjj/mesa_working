      subroutine laginf (nlag,lstlag,nsts,nflag)
      implicit integer(a-z)
      character *4 itoc
      dimension nlag(nsts), lstlag(nflag,nsts)
      call iosys ('read integer "no. lag orbs" from rwf',nsts,nlag,0,0)
      do 10 i=1,nsts
      call iosys ('read integer "lag orbs'//itoc(i)//'" from rwf',nlag(i
     1 ),lstlag(1,i),0,0)
   10 continue
      return
      end
