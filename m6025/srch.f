*deck srch
      subroutine srch(unit6,unit9,debug)
      implicit real*8 ( a - h, o - z )
      integer unit6,unit9
      common/mprin/l
      common/irg/irang(30)
      common/kaid/aid(14,14),aid2(14,14),naid
c  naid must be less than or equal to numax
      common/kar/u(30),usav(30),grad(30),z(30),zsav(30),bl(30),bu(30),v(
     130),vsav(30),gradu(30),nuflag(30),nbd(30)
      common /kchg/ gsav( 300,14),b( 300),bsav( 300),nersdim
      dimension gg( 300,14)
      common/kvar/nu,nprob,ners,numax,nerspnu,phi,solves,ngrad,
     1netasr,gradsq,zsq
      dimension eta3(3),phi3(3),etas2(50)
      common/d1/ yz,timesum,nbdsum,nuvar  ,rus ,iter
      common/doyle/nprobs,iters,grdzmin,phimin,relphi,zsqmin,np0,nprint1
     1,nprint2,netasch,mono,ntable,npiv,apzero,cmin,cmax,relpphi,reldgmx
     1,etavrel,etahalt,usqmx ,mxfcn,icrow
     1,dmult
      logical solves, debug
c
c
      save dgmxpre,emult,exp1
c
      common /pb/ p
      data naid/14/
      data lshift/5/
      if(netasr-2)116,110,118
c     netasr = 2, read data
  110 continue
      cmaxs=cmax
      cmins=cmin
      icrset=0
      emult=dmult
      net=netasch
      if (nuvar .eq. 1)  icrow = 0
      netasr=1
      sum =2*ntable
      exp1= sum-2.0d0
      exp1=1.0d0/exp1
      tlstsq=0.0d0
      dgmxpre=1.0d+99
  116 ngrad=1
      nbr1=1
  118 storen=1.0d0
      netasch=net
      igls=0
      nkp=0
      do 11190 k=1,nuvar
11190 usav(k)=u(k)
11191 if(nbdsum)11192,11192,11196
c     nbdsum=0, no bds
11192 ind=0
      nbdsum=0
      do 11194 k=1,nu
      if(nuflag(k)) 11194,11193,11194
11193 ind=ind+1
      v(k)=u(ind)
11194 continue
      call phig
      go to 11197
11196 call bnds(unit6,unit9,debug)
      if(nuvar-1) 11680,1000,11197
 1000 np0=0
      igls=0
      icrow = 0
      mono = 1
      ntable = 2
11197 nerpvar=ners+nuvar
      if(nersdim-nerpvar)11198,11199,11199
11198 if (debug) write(l,59)nerpvar,nersdim
      go to 11313
11199 nersp1=ners+1
      go to (11200,11250,11208,11470) nbr1
11200 gradsq=0.0d0
      som=0.0d0
      do 11207 k=1,nuvar
      b(k)=0.0d0
      do 11201 ner=1,ners
11201 b(k)=b(k)+gsav(ner,k)*gsav(ner,k)
      if(b(k)-som) 11217,11217,11203
11203 som=b(k)
11217 if(nbdsum) 11209,11215,11209
11209 gradsq=gradsq+gradu(k)*gradu(k)
      go to 11207
11215 gradsq=gradsq+grad(k)*grad(k)
11207 continue
11208 phistrt = phi
      fipprev = phi
11210 diagmax=som
      diag1=cmin*som
      diag2=cmax*som
      pfac = cmax/(cmin+1.d-8)
      pfac=pfac**exp1
      pfacsq=pfac**2
      pn1 =sqrt(diag1)/pfac
c     write(unit9,7)som,pfac
c     write(unit9,49)(b(k),k=1,nuvar)
      fimin=1.0d+250
      if(np0)11211,11211,11212
11211 npstep=0
      p=pn1
      go to 11213
11212 npstep=-1
      p=0.0d0
11213 ngrad=0
      nbr1=2
11214 p=p*pfac
      npstep=npstep+1
      if(igls) 10204,10204,10205
10204 igls=1
      do 11227 ner=1,ners
      b(ner)=bsav(ner)
      do 11227 k=1,nuvar
11227 gg(ner,k)=gsav(ner,k)
      if(icrow) 122,122,120
  120 call normal(gsav,nersdim,ners,nuvar,z,1)
      call rownorm(gsav,nersdim,b,ners,nuvar)
  122 continue
      do 126 i=1,nuvar
      do 124 k=1,nuvar
  124 aid(i,k)=0.d0
  126 aid(i,i)=1.d0
      call glss(ners,nuvar,1,nuvar,nrank,gsav,nersdim,b,nersdim,aid,naid
     1,z,nu)
      do 127 i=1,nuvar
      b(i)=sdot(ners,gsav(1,i),1,bsav(1),1)
  127 bsav(i)=b(i)
      call sgefa(aid,naid,nuvar,gsav(1,1),info)
      call sgedi(aid,naid,nuvar,gsav(1,1),gsav(1,2),gsav(1,3),0)
10205 if(p-1.d-100) 11225,11225,11216
11216 do 11219 i=1,nuvar
      b(i)=bsav(i)
      do 11218 k=1,nuvar
      aid2(i,k) = 0.d0
11218 gsav(i,k)=aid(i,k)
11219 aid2(i,i)=1.d0
      n1=nuvar+1
      n2=2*nuvar
      do 11222 i=n1,n2
      b(i)=0.d0
      do 11222 k=1,nuvar
      gsav(i,k)=0.d0
      if(icrow) 11223,11223,11224
11223 gsav(k+nuvar,k)=p
      go to 11222
11224 gsav(k+nuvar,k)=expad(p,-irang(k))
11222 continue
      call glss(n2,nuvar,1,nuvar,nrank,gsav,nersdim,b,nersdim,aid2,naid,
     1z,nu)
11225 ndefic=nuvar-nrank
      if(icrow) 132,132,130
  130 call normal(gsav,nersdim,ners,nuvar,z,2)
  132 continue
      if(ndefic)11228,11230,11228
11228 write(unit9,38)npstep,ndefic
      if(nuvar-1) 11680,11290,11290
11230 zsq=0.0d0
      do 11231 k=1,nuvar
11231 zsq=zsq+z(k)*z(k)
      if(nbdsum) 11232,11233,11232
11232 yz=sdot(nuvar,gradu,1,z,1)
      go to 11234
11233 yz=sdot(nuvar,grad,1,z,1)
11234 if(abs(yz)-grdzmin)11237,11237,11235
11235 rtgrad=sqrt(gradsq)
      rtzsq=sqrt(zsq)
      sum=rtzsq*rtgrad
      costhet=-yz/sum
11237 do 11245 k=1,nuvar
11245 u(k)=usav(k)+z(k)
      go to 11191
11250 sum=(phistrt-phi)/phistrt
      rus=sum
      if(nprint1) 11260,11260,11255
11255 if (debug) write(l,11)npstep,p,rtzsq,costhet,phi,sum,tlstsq
11260 if(nuvar-1) 11266,11261,11261
11261 if(phi .ge. phistrt .or. phi .ge. fimin) go to 11290
11256 if( abs(sum) .ge. 1.0d-25) go to 11265
      if(icrset) 11305,11305,11257
11257 write(unit9,60)
      nuvar=-nuvar
      go to 11650
11265 if(p .ne. 0.d0 .or. iter .ge. 3 ) go to 11266
      do 11267 k=1,nu
      if(nbd(k)) 11267,11267,11268
11268 if(nbd(k)-2) 11273,11271,11271
11271 if(bu(k)-v(k)-1.0d-4) 11272,11272,11273
11272 phi=1.0d100
      go to 11290
11273 if(nbd(k)-2) 11274,11274,11267
11274 if(v(k)-1.0d-4-bl(k)) 11272,11272,11267
11267 continue
11266 fimin=phi
      zmin=rtzsq
      zsqopt=zsq
      yzopt=yz
      popt=p
      nmin=npstep
      do 11270 k=1,nuvar
11270 zsav(k)=z(k)
      do 11275 k=1,nu
11275 vsav(k)=v(k)
      if(nuvar-1) 11313,11276,11276
11276 if(sum-relpphi)11290,11290,11285
11285 sum1=cmin
      sum2=cmax
      if(npstep-ntable) 11282,11282,11278
11278 sum2=cmax*pfacsq
11282 if(npstep-1) 11286,11286,11284
11284 sum1=cmin*pfacsq
11286 continue
      cmin=sum1
      cmax=sum2
      netasch=0
      go to 11316
11290 if(p)11292,11292,11293
11292 p=pn1
11293 if(mono)11302,11302,11294
11294 if(phi-fipprev)11296,11298,11298
11296 storen=-1.0d0
      go to 11302
11298 storep=1.0d0
      sum3=storen*storep
      if(sum3)11300,11302,11302
11300 if(fimin-phistrt)11307,11302,11302
11302 fipprev=phi
      if(npstep-ntable)11214,11303,11303
11303 if(mono)11316,11316,11304
11304 if(fimin-phistrt)11306,11305,11305
11305 antable=ntable
      if(icrset .eq. 1.or. nkp .lt. 2) go to 10309
      cmax=cmaxs
      cmin=cmins
      nkp=0
10306 icrset=1
      if(icrow) 10307,10307,10308
10307 icrow=1
      go to 10310
10308 icrow=0
10310 if (debug) write(l,64) icrow
      igls=0
      do 1038 i=1,ners
      do 1038 k=1,nuvar
 1038 gsav(i,k)=gg(i,k)
      go to 11210
10309 continue
      if(nkp-3) 11301,11301,11311
11311 write(unit9,100)
  100 format(' stop 11305 ', 50x,'stop 11305 ',45x,'stop 11305')
11313 iter=99999
      do 11315 k=1,nuvar
11315 u(k)=usav(k)
      go to 11680
11301 sum1=cmax*pfacsq
      sum=pfacsq**ntable
      sum=cmax*sum
c     write(unit9,33)cmin,sum1,cmax,sum
      cmin=sum1
      cmax=sum
      npstep=0
      nkp=nkp+1
      go to 11214
11306 sum=cmax*pfac*pfac
c     write(unit9,29)cmax,sum
      cmax=sum
      go to 11312
11307 if(npstep-ntable)11308,11309,11309
11308 sum=cmax/pfacsq
c     write(unit9,29)cmax,sum
      cmax=sum
11309 if(npstep-3)11310,11316,11312
11310 sum=cmin/(pfac*pfac)
      go to 11314
11312 sum=cmin*pfacsq
11314 cmin = sum
11316 do 11320 k=1,nu
11320 v(k)=vsav(k)
      do 11318 k=1,nuvar
11318 u(k)=usav(k)+zsav(k)
      sum1=popt*popt
      sum2=pfac*pfac
11332 if (debug) write(l,14)popt,phistrt,fimin
c11332 write(unit9,13)nmin,popt,phistrt,fimin,zmin,yzopt
c     write(unit9,49)(zsav(k),k=1,nuvar)
c     write(unit9,15)gradsq,yz,timesum
      g1=0.d0
      do 11333 k=1,nuvar
11333 g1=g1+grad(k)**2
      if (debug) write(l,16) g1,(grad(k),k=1,nuvar)
c     write(unit9,49)(grad(k),k=1,nuvar)
      if(fimin-phistrt)11336,11334,11334
11334 write(unit9,43)fimin,phistrt
      stop
11336 phi=fimin
      zsq=zsqopt
      yz=yzopt
      if(phi-1.d-24) 11680,11680,11340
11340 if(netasch)11650,11650,11430
c     netasch=1, make eta search
11430 phi3(1)=phistrt
      eta3(1)=0.d0
      phi3(2)=fimin
      eta3(2)=1.d0
      eta=2.d0
      eta3(3)=eta
      nbr2=1
      nbr1=4
      nfeval=0
      nrt=-1
      nvx=0
      netas2=0
11440 nfeval=nfeval+1
      netas2=netas2+1
      etas2(netas2)=eta3(2)
      if(nfeval-mxfcn)11445,11445,11640
11445 do 11450 k=1,nuvar
11450 u(k)=usav(k)+eta*zsav(k)
      go to (11191,11644) nbr2
11470 if(phi-phi3(2))11480,11570,11570
11480 if(eta-eta3(2))11484,11484,11482
11482 phi3(1)=phi3(2)
      eta3(1)=eta3(2)
      go to 11486
11484 phi3(3)=phi3(2)
      eta3(3)=eta3(2)
11486 phi3(2)=phi
      eta3(2)=eta
c has monotonicity been broken (n0 is lt 0)
      if(nrt)11488,11640,11582
11488 eta=emult*(eta+1.d0)
      nrt=nrt-1
      if(lshift+nrt) 11489,11440,11440
11489 eta=eta/emult-1.d0
      emult=emult*2.d0
      go to 11640
11570 nrt=1
      emult=dmult
11571 if(nvx)11572,11574,11576
11574 if(eta-eta3(2))11580,11580,11578
11578 phi3(3)=phi
      eta3(3)=eta
      go to 11582
11580 phi3(1)=phi
      eta3(1)=eta
11582 continue
11581 etanum=0.5d0*(phi3(1)*(eta3(2)**2-eta3(3)**2)-phi3(2)*(eta3(1)**2-
     1eta3(3)**2)+phi3(3)*(eta3(1)**2-eta3(2)**2))
      etaden= phi3(1)*(eta3(2)-eta3(3))-phi3(2)*(eta3(1)-eta3(3))+phi3(3
     1)*(eta3(1)-eta3(2))
      if(abs(etaden)-1.0d-30)11640,11640,11583
11583 etavrtx=etanum/etaden
      eta=etavrtx
      if(eta-eta3(2))11584,11640,11586
11584 nvx=-1
      go to 11590
11586 nvx=+1
      go to 11590
11572 eta=.5d0*(eta3(3)+eta3(2))
      phi3(1)=phi
      eta3(1)=etavrtx
      nvx=0
      go to 11590
11576 eta=0.5d0*(eta3(2)+eta3(1))
      phi3(3)=phi
      eta3(3)=etavrtx
      nvx=0
11590 releta=abs((eta-eta3(2))/eta3(2))
      if(nprint2)11532,11532,11530
11530 write(unit9,23)
      write(unit9,35)(eta3(k),k=1,3)
      write(unit9,24)
      write(unit9,35)(phi3(k),k=1,3)
11532 if(releta-etavrel)11640,11640,11534
11534 if(abs(phi3(1)+phi3(3)-2d0*phi3(2))/
     >        (phi3(1)+phi3(2)+phi3(3))-.05d0) 11640,11440,11440
11640 phi=phi3(2)
      etavrtx=eta3(2)
      eta=eta3(2)
      nbr2=2
      go to 11445
11644 sum=(phistrt-fimin)/phistrt
      sum1=etahalt*sum
      sum=(fimin-phi)/fimin
      if (debug) write(l,25)fimin,phi,etavrtx,sum
c     write(unit9,44) nfeval
   44 format(20x,'no of function evaluations in eta search=',i4)
c     write(unit9,35)(etas2(k),k=1,netas2)
      if(sum)11641,11642,11642
11641 if (debug) write(l,36)phi,fimin
      phi=fimin
      eta=1.0d0
11642 if(sum-sum1)11643,11643,11650
11643 if(sum.lt.1.d-200.or.iter.le.3.or.nfeval.lt.3)go to 11650
      netasch=0
      if (debug) write(l,32)sum,sum1
11650 sum=abs((diagmax-dgmxpre)/dgmxpre)
      if(sum .gt. reldgmx .or. iter .lt. 7 ) go to 11670
11660 ngrad=-1
      netasr=3
      if (debug) write(l,22)dgmxpre,diagmax,sum,reldgmx
      nbr1=3
      gradsq=1.0d+99
      if (debug) write(l,37)
11670 continue
      dgmxpre=diagmax
11680 return
    1 format(i4)
    2 format('    srch 11320, unconstrained par vec u =')
    3 format(e8.2)
    4 format(' nprint1 =',i4,', 0(1) to by-pass(print) p-table at 11255'
     c)
    6 format(' relpphi =',1pe9.2, ', exit from p-srch when rel phi drop'
     c' greater than relpphi')
    7 format('    srch 11207, diagmax =',1pe13.5,
     c', pfac =',1pe9.2,', diagonal terms are')
    8 format('    srch 11194, phi =',1pe22.14,', nbdsum =',i3,', par vec
     c v(k) =')
   11 format(' n=',i3,', p=',1pe9.2,', z=',1pe9.2,', cos=',1pe9.2,', phi
     c=',1pe21.14,', rel phi drop =',1pe9.2                 ,/' lstsq ti
     cme =',1pe11.4)
   12 format('    srch 11285, rel phi drop =',1pe10.2,', is grtr than re
     clpphi =',1pe10.2,   /'psrch exit, cmax= ',1pe15.8,'now=',1pe15.8,'
     1cmin=',1pe15.8,'now=',1pe15.8)
   13 format('    srch 11332, opt n =',i3,', opt p=',1pe12.5,', start ph
     ci =',1pe22.14,', end phi =',1pe22.14,/'  opt z =',1pe9.2,', grad.z
     c optimal =',1pe9.2,                                           ', p
     cfac  sq =',1pe9.2,', popt sq =',1pe9.2,', z vec =')
   14 format('      srch11332,opt p=',1pe12.5,',start phi=',1pe22.14,',e
     1nd phi=',1pe22.14)
   15 format('    srch 11332, gradsq =',1pe10.2,', z.grad =',1pe10.2,','
     c' phig    timesum =',1pe9.2,', gradient =')
   16 format('   gradsq=',1pe10.2,'gradient vector=',1p5e13.5,/(2x,1p7e1
     13.5))
   17 format(' mono    =',i4,', 0(1) to continue(halt) table at min')
   18 format(' cmin    =',1pe9.2,', starting p sq = cmin(diagmax)')
   19 format(' cmax    =',1pe9.2,',      end p sq = cmax(diagmax)')
   20 format(' reldgmx =',1pe9.2,', halt differentiation and p-table whe
     1n relative change in diag  max is less than reldgmx')
   21 format(' netasch =',i4,', 0(1) to omit(include) opt z(p) eta searc
     ch')
   22 format('    srch 11660, previous max diag term =',1pe12.5,' and pr
     1esent diag term =',1pe12.5,' have rel change =',1pe10.3,/'  less t
     2han input reldgmx=',1pe9.2,',halt differentiation')
   23 format('    srch 11530, parabola min  sequence= '  )
   24 format('    srch 11530, phi triplet =')
   25 format('    srch 11640, p-srch opt phi =',1pe22.14,', eta-srch opt
     1 phi =',1pe22.14,', opt eta =',1pe11.3,/' relative phi drop =',1pe
     210.2)
   26 format(' etavrel =',1pe9.2,', eta-srch exit when par vrtx rel ch l
     cess th etavrel')
   27 format(' nprint2 ',i4,',0(1) to by-pass(print)    srch triplets')
   28 format(' np0     =',i4,', 0(1) to by-pass(include) p=0 in p-tabl')
   29 format('    srch 11307, old cmax =',1pe10.3,' now reset to ',1pe10
     c.3)
   30 format('    srch 11314, old cmin =',1pe10.3,' now reset to ',1pe10
     c.3)
   31 format(' etahalt =',1pe9.2,', halt eta-srch when rel phi drop less
     c  than etahalt(rel p-table phi drop)')
   32 format('    srch 11642, rel eta-srch phi drop =',1pe11.3,' is less
     c than etahalt(p-table rel phi drop)=',1pe11.3,', halt eta-srch')
   33 format('    srch 11305, continue p-table with old cmin=',1pe10.3,'
     c reset to ',1pe10.3,', old cmax=',1pe10.3,' reset to',1pe10.3)
   35 format(1p5e22.14)
   36 format('    srch 11641,    srch phi =',1pe21.14,' is greater than'
     c' p-srch opt phi =',1pe21.14,', p-srch opt phi restored')
   37 format('011660, differentiation halted, symbol gradsq=1.0d+99')
   38 format('    srch 11228, singular system, npstep =',i3,', ndefic ='
     1,i3,', pivot1(k) =')
   39 format('    srch 11228, pivot2(k) =')
   40 format('    srch 11228, jcol(k) =')
   41 format(30i4)
   42 format('    srch 11229, error halt, sing syst for const p =',1pe10
     c.3)
   43 format('    srch 11334, error halt, end phi =',1pe22.14,' is great
     cer than start phi =',1pe22.14)
   49 format(1p10e13.5)
   56 format(' apzero  =',1pe9.2,', approximate 0 for lstsq rank determi
     cnation')
   58 format(' ntable  =',i4,', number of p-table entries')
   59 format('    srch11198,nerpvar =',i6,' exceeds control input ners '
     c' dim =',i6,', make nersdim = nerpvar in control order and dimensi
     cons b, g')
   60 format('    srsh  11256rel phi drop in p search =0 exit check conv
     1.  ')
   62 format(' icrow= ',i5,3x,'1(0) col followed by row max scaling (not
     1)  desired  ')
   64 format(' scaling parameter icrow has  been set to',i5,3x,'and p ta
     1ble restarted ')
   65 format(' mxfcn= ',i5,3x,'max no. of fcn eval-1   allowed in etasch
     '')
   66 format(' dmult= ',f10.4,3x,' eta multiplyer in etasrch')
   74 format('    srch 11196, ner =',i4,',           nuvar = ',i3,', eq'
     c' =',1pe10.3,', eqj(k) =')
   96 format(1p10e13.5)
  108 format(' npiv    =',i4,', 0(1) when max modified col vec (modified
     c length/orig length) determines next pivot')
c     subrt srch ends here
      end
