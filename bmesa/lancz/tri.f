*deck Program Tri.f 
c***begin prologue     tri
c***date written       960718   (yymmdd)
c***revision date               (yymmdd)
c***keywords           lanczos
c***author             schneider, b. i.(nsf)
c***source             tri
c***purpose            
c***                   
c***references
c***routines called    iosys, util and mdutil
c***end prologue       tri
c
      program tri
c
      implicit integer (a-z)
      parameter ( maxreg=100, mgr=1 )
      character*4096 ops
      character*2 itoc
      character*80 cpass, chrkey
      character*24 coord, key, typpot
      character*800 card
      character*128 filbec
      logical prn, dollar, logkey, reuse
      real*8 edge, lanc, ham 
      integer*8 phamil
      dimension edge(maxreg+1)
      dimension  prn(40)
      common/io/inp, iout      
      pointer (pl,lanc(1))
      pointer (phamil,ham(1))
c
      call drum
      write(iout,*)
      call iosys ('read character options from rwf',-1,0,0,ops)
      dim=1
      write(iout,1)
c
      call iosys ('read character "bec filename" from rwf',
     1             -1,0,0,filbec)
      call iosys ('open bec as new',0,0,0,filbec)
      call iosys('rewind all on bec read-and-write',0,0,0,' ')
      coord=chrkey(ops,'dimension-1','x',' ')
      iter=intkey(ops,'number-of-lanczos-iterations',10,' ')
      ngrds=1
      do 10 i=1,20
         prn(1)=.false.
 10   continue   
      call basis(ops,cpass,card,i,phamil,edge,ngrds,coord,
     1              key,typpot,n,nphy,nwham,reuse,prn)
      v=1
      a=v+(iter+1)*nphy
      b=a+iter
      scr=b+iter
      vec=scr+nphy
      need=wpadti(vec+iter*iter)
      call memory(need,pl,ngot,'lancz',0)
      call initv(lanc(v),nphy)
      call lancz(lanc(v),ham,lanc(a),lanc(b),lanc(vec),
     1           lanc(scr),nphy,iter)
      call memory(-ngot,pl,idum,'lancz',idum)
      call chainx(0)               
      stop
 1    format(/,20x,'lanczos tridiagonalization')      
      end









