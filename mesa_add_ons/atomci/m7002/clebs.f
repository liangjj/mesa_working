*deck @(#)clebs.f	2.1  10/10/91
      program clebs
      implicit integer(a-z)
      common // z(1)
      real*8 z
      logical logkey, prnt
      dimension a(1)
      character *4096 ops
      character *128  filci
      character *3 ans
      equivalence (a(1),z(1))
      common /memory/ ioff
      common /io/ inp, iout
c
c
      call drum
      call iosys ('read character options from rwf',-1,0,0,ops)
      call iosys ('read character "atomic ci filename" from rwf',-1,
     1             0,0,filci)
      call iosys ('open atomci as unknown',0,0,0,filci)
      lmax=intkey(ops,'maximum-l-value',8,' ')
      prnt=logkey(ops,'print=m7002',.false.,' ')
      bufsiz=intkey(ops,'buffer-size',100000,' ')
      l1=lmax+1
      l2=2*lmax+1
      need=l1*l1*l1*l2*l2
      bufsiz=min(bufsiz,need)
      need=need+iadtwp(need)
      minwds=bufsiz+iadtwp(bufsiz)
      words=min(need,minwds)+1000 
      call iosys ('read integer maxsiz from rwf',1,maxcor,0,' ')
      if (words.gt.maxcor) then
          call lnkerr('not enough memory for calculation:quit')
      endif
      call iosys ('write integer maxsiz to rwf',1,words,0,' ')
      call getscm(words,z,ngot,'m7002',0)
      cleb=ioff
      acleb=wpadti(ioff)
      icleb=wpadti(cleb+bufsiz)
      call iosys('write integer "cg buffer size" to atomci',1,
     1            bufsiz,0,' ')
      call iosys('does "cg coef" exist on atomci',0,0,0,ans)
      if (ans.ne.'yes') then
          call iosys('create integer "cg coef" on atomci',-1,0,0,' ')
      endif
      call racah(z(cleb),a(acleb),a(icleb),lmax,bufsiz,nowrts,last)
      if(ans.ne.'yes') then
         call iosys('endfile "cg coef" on atomci',0,0,0,' ')
      endif
      call iosys('write integer "no. cg writes" to atomci',1,nowrts,
     1            0,' ')
      call iosys('write integer "size of last cg write" to atomci',1,
     1            last,0,' ')
      call iosys('rewind all on atomci read-and-write',0,0,0,' ')
      if (prnt) then
          call rdcleb(z(cleb),a(acleb),a(icleb))
          call iosys('rewind all on atomci read-and-write',0,0,0,' ')
      endif
      call iosys('close atomci',0,0,0,' ')
      call iosys ('write integer maxsiz to rwf',1,maxcor,0,' ')
c
c
      call chainx(0)
      stop
      end




