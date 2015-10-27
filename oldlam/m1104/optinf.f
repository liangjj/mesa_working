      subroutine optinf (nopt,lstopt,kepvnl,kpdiag,nsts,nfopt)
      implicit integer(a-z)
      character *4 itoc
      dimension nopt(nsts), lstopt(nfopt,nsts)
      dimension kpdiag(nsts,nsts)
      dimension kepvnl(nfopt,nsts)
      call iosys ('read integer "no. opt orbs" from rwf',nsts,nopt,0,0)
      do 10 i=1,nsts
      call iosys ('read integer "kpdiag'//itoc(i)//'" from rwf',nsts
     1 ,kpdiag(1,i),0,0)
      call iosys ('read integer "opt orbs'//itoc(i)//'" from rwf',nopt(i
     1 ),lstopt(1,i),0,0)
      call iosys ('read integer "kepvnl'//itoc(i)//'" from rwf',nopt(i)
     1 ,kepvnl(1,i),0,0)
   10 continue
      return
      end
