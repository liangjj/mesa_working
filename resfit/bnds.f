      subroutine bnds
      common/kar/u(30),usav(30),grad(30),z(30),zsav(30),bl(30),bu(30),v(
     130),vsav(30),gradu(30),nuflag(30),nbd(30)
      common /kchg/ gsav( 300,14),b( 300),bsav( 300),nersdim
      common/kvar/nu,nprob,ners,numax,nerspnu,phi,solves,ngrad,
     1netasr,gradsq,zsq
      common/d1/ yz,timesum,nbdsum,nindep1,rus,iter
      common/doyle/nprobs,iters,grdzmin,phimin,relphi,zsqmin,np0,nprint1
     1,nprint2,netasch,mono,ntable,npiv,apzero,cmin,cmax,relpphi,reldgmx
     1,etavrel,etahalt,usqmx ,mxfcn,icrow
     1,dmult
      common /pb/ p
      data eps/1.e-6/
      irdb = 0
      if(ngrad-2)88,34,48
34    write(9,3)
      lter = 0
      write(9,1)(nbd(k),k=1,nu)
      usqd=1./usqmx
      rus=-1.
      nbdsum=0
      do 36 k=1,nu
      if(nuflag(k))36,36,35
c     nuflag(k)=1, v(k) fixed
   35 nbd(k)=0
   36 nbdsum=nbdsum+nbd(k)
      write(9,4)nbdsum
      write(9,1) (nbd(k),k=1,nu)
      if(nbdsum)37,37,39
c  nbdsum =0,  no bounds
   37 ind=0
      do 38 k=1,nu
      if(nuflag(k)) 41,40,41
   40 ind=ind+1
      u(ind)=v(k)
   41 continue
   38 continue
   39 return
c     compute unbounded u(k) to yield starting bounded v(k)
c     nbdsum non 0, parameters bounded, nbd(k) = (0, no bds), (1, l bd o
c     nly), (2, l and u bds), (3, u bd only)
   48 ind=0
      do 80 k=1,nu
      if(nuflag(k)) 80,49,80
   49 ind=ind+1
      if(nbd(k))52,50,52
c     nbd(k)=0, v(k) unbounded or fixed
   50 u(ind)=v(k)
      go to 80
c     nbd(k) non 0
   52 if(nbd(k)-2)54,60,74
c     nbd(k)=1, l bd only
   54 vma=v(k)-bl(k)
      if(vma)56,56,58
c     v(k) below or on l bd
   56 write(9,5) k,v(k),bl(k)
      v(k) = v(k) + eps
      vma  = vma + eps
      if (vma .lt. 0.)   call exit
c     nbd(k)=1, l bd only, v(k) above l bd
   58 u(ind)=1.0/sqrt(vma)
      go to 80
c     nbd(k)=2, l and u bds
   60 vma=v(k)-bl(k)
      sum1=0.5*(bu(k)+bl(k))
      if(vma)62,62,64
c     l and u bds, v(k) below or on l bd
   62 write(9,10) k,v(k),bl(k),bu(k)
c     l and u bds, v(k) above l bd
      v(k) = v(k) + eps
      vma = vma + eps
      if (vma .lt. 0.)   call exit
   64 vmb=v(k)-bu(k)
c     l and u bds, v(k) between
      if (vmb .lt. 0.)   go to 68
      write (9,10) k,v(k),bl(k),bu(k)
      v(k) = v(k)-eps
      vmb  = vmb -eps
      if (vmb .gt. 0.)   call exit
   68 sum2=(bu(k)-bl(k))*.5
      sum1=v(k)-sum1
      sum2=sum1/sum2
      u(ind)=gtan(sum2)
      go to 80
c     nbd(k)=3, u bd only
   74 bmv=bu(k)-v(k)
      if(bmv)76,76,78
c     u bd only, v(k) above or on u bd
   76 write(9,11) k,vk,bu(k)
c     u bd only, v(k) below u bd
      bmv = bmv + eps
      v(k) = v(k) - eps
      if (bmv .lt. 0.)   call exit
   78 u(ind)=1.0/sqrt(bmv)
   80 continue
      if (irdb.eq.1)  go to 160
      write(6,7)
      write(6,8) (u(k),k=1,nindep1)
      return
c     nbdsum non 0, param bounds
   88 ind=0
      do 112 k=1,nu
      if(nuflag(k)) 112,89,112
   89 ind=ind+1
      if(nbd(k))92,90,92
c     nbd(k)=0, v(k) unbounded or fixed
   90 v(k)=u(ind)
      go to 112
c     nbd(k) non 0, v(k) bounded
   92 uksq=u(ind)**2
      if(nbd(k)-2)94,106,110
c     nbd(k)=1, l bd only
   94 v(k)=bl(k)+1.0/uksq
      go to 112
c     nbd(k)=2, l and u bds
  106 sum=u(ind)
      v(k)=0.5*(bu(k)+bl(k) +(bu(k)-bl(k))*gatan(sum))
      go to 112
c     nbd(k)=3, u bd only
  110 v(k)=bu(k)-1.0/uksq
  112 continue
      if (ngrad .ne. 1 )           go to 160
  167 ind=0
      do 160 k=1,nu
      if(nuflag(k)) 160,168,160
  168 ind=ind+1
      if(nbd(k)-2) 170,169,170
  169 us=1.-abs(gatan(u(ind)))
      if(us-usqd) 171,171,160
  171 if(rus .lt. 0.)  go to 160
      lter=lter+1
      if ( lter .lt. 13 )          go to 160
      nindep1 = nindep1 - 1
      lter = 0
   16 format(* bnds 171 v(*,i2,*) = *,e20.12,*nbdsum=*,e20.12,*us=*,e20.
     112,*no.  indep var red. to *,i3  )
      write(9,16) k,v(k),nbdsum,us,nindep1
      go to 172
  170 uksq=u(ind)**2
      lter = lter + 1
      if ( lter .lt. 13 )          go to 160
      lter = 0
      if(uksq .lt. usqmx .or. rus .lt. 0. ) go to 160
      nindep1=nindep1-1
      write(9,15) k,v(k),nbdsum,uksq,nindep1
   15 format(* bnds 160 v(*,i2,*) =*,e20.12,*nbdsum=*,e20.12,*uksq=*,e20
     1.12,*no. indep var red. to *,i3 )
  172  nbdsum=nbdsum-nbd(k)
       nuflag(k)=1
       nbd(k)=0
      irdb = 1
      go to 48
  160  continue
      if (irdb.eq.0) go to 155
      do 151  k=1,nindep1
  151 usav(k) = u(k)
  155 continue
      if(nbdsum .eq. 0) nbdsum=-1
      call phig
      if(ngrad)154,154,120
  120 ind=0
      sum1=0.0
      do 150 k=1,nu
      if(nuflag(k)) 150,122,150
c     nuflag(k)=0, v(k) varied
  122 ind=ind+1
       if(nbd(k))163,164,166
  163  call exit
  164  sum=1.
       go to 130
  166 sum=u(ind  )**2
      nbdkm2=nbd(k)-2
      if(nbdkm2)124,128,124
c     l or u bnd only
  124 sum=u(ind  )*sum
      sum=2.0/sum
      if(nbdkm2)126,128,130
c     l bd only
  126 sum=-sum
      go to 130
c     nbd(k)=2, l and u bds
  128 sum=(1.0+abs(u(ind  )))**2
      sum=.5*(bu(k)-bl(k))/sum
  130 gradu(ind)=grad(ind)*sum
      sum1=sum1+gradu(ind)**2
      do 140 ner=1,ners
  140 gsav(ner,ind)=gsav(ner,ind)*sum
  150 continue
      write(6,12)sum1
      write(6,14)(gradu(k),k=1,nindep1)
  154 return
    1 format(36i2)
    2 format(e8.2)
    3 format(* input parameter bounds flag vector nbdk(k) =*)
    4 format(* bnds 36, nbdsum =*,i3,*, par bd vec, after possible modif
     cication by par flag vec, =*)
    5 format(* bnd 56, for param index k =*,i3,* v(k) =*,1pe22.14,* is l
     c th or = l bd bl(k) =*,1pe14.6 )
    7 format(* bnds 80, u(k) to produce starting v(k) is*)
   11 format(* bnds 76, u bd only, v(k) above u bd, k =*,i3,*, v(k) =*,1
     cpe22.14,*, bu(k) =*,1pe14.6  )
    8 format(1p8e16.6)
   10 format(* bnds 62, l and u bnds, v(k) beyond bnds, k =*,i3,*, v(k)*
     c* =*,1pe22.14,*, bl(k) =*,1pe14.6,/*  bu(k) =*,1pe14.6  )
   12 format(* bounds 150, gradusq =*,1pe10.3,*, u-gradient vec =*)
   14 format(1p10e12.3)
c     subrt bnds ends here
      end
