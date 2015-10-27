*deck opac
      subroutine opac
      implicit real *8 (a-h,o-z)
      common/mprin/l
      common/kaid/aid(14,14),aid2(14,14),naid
      common/kar/u(30),usav(30),grad(30),z(30),zsav(30),bl(30),bu(30),v(
     130),vsav(30),gradu(30),nuflag(30),nbd(30)
      common /kchg/ gsav( 300,14),b( 300),bsav( 300),nersdim
      common/kvar/nu,nprob,ners,numax,nerspnu,phi,solves,ngrad,
     1netasr,gradsq,zsq
      logical solves
      common/doyle/nprobs,iters,grdzmin,phimin,relphi,zsqmin,np0,nprint1
     1,nprint2,netasch,mono,ntable,npiv,apzero,cmin,cmax,relpphi,reldgmx
     1,etavrel,etahalt,usqmx ,mxfcn,icrow
     1,dmult
      common/d1/ yz,timesum,nbdsum,nindep1,rus,iter
      data dmult/1.d0/
      data usqmx/1.d06/
      data mxfcn/3/
      data icrows/1/
      data nprobs/1/
      data iters/10 /
      data grdzmin/1.d-20/
      data phimin/1.d-25/
      data relphi/1.d-13/
      data zsqmin/1.d-24/
      data np0s/1/
      data nprint1/0/
      data nprint2/0/
      data net/1/
      data mono/1/
      data ntable/3/
      data npiv/0 /
      data apzero/1.d-14/
      data cmins/1.d-3/
      data cmaxs/1.d-0/
      data relpphi/.8d0/
      data reldgmx/1.d-50/
      data etavrel/1.5d-02/
      data etahalt/1.0d-02/
      data kread/1/
      data nersdim/300/
      solves=.false.
  170 ngrad=2
  480 nindep1=0
      do 460 k=1,nu
      if(nuflag(k)) 460,470,460
c           nuflag (k)=0, vary u(k)
  470 nindep1=nindep1+1
  460 continue
      netasch=net
      cmin=cmins
      cmax=cmaxs
      np0=np0s
      rus=-1.d0
      icrow=icrows
      call phig
      call bnds
      write(9,5)iters
        write(9,12)relphi
      write(9,11)phimin
      write(9,61)grdzmin
      write( 9,18)zsqmin
  180 call second(tstart)
      write(9,51)nprob
      write(9,55)
      write(9,56)
      write(9,57)
      write(9,58)
c     print starting par vec v(k), flags, bounds
      write(9,52)
      do 182 k=1,nu
  182 write(9,53)nuflag(k),nbd(k),bl(k),v(k),bu(k),k
c     compute unbounded u(k) to yield initial v(k)
      ngrad=3
      call bnds
      if(nindep1) 235,235,183
  183 iter=0
      l=9
      phi=0.0d0
  184 netasr=2
      call srch
      write(l,100) (v(i),i=1,nu)
  186 phiprev=1.0d+99
      netasr=1
      l=6
  195 iter=iter+1
      go to 240
  230 solves =.true.
  231 write(9,4)
      if(nindep1) 232,233,233
  232 nindep1=-nindep1
  233 continue
      ngrad=1
      if (nbdsum) 254,254,256
  254 ind = 0
      do  258  k=1,nu
      if (nuflag(k)) 258,257,258
  257 ind = ind+1
      v(k) = u(ind)
  258 continue
      call phig
      go to 235
  256 call bnds
  235 write(9,8)phi
      write(9,105) iter
c     statistical estimator normalized when practical
      do 620 i=1,nindep1
      do 610 k=1,nindep1
      aid(i,k)=0.d0
      do 600 j=1,ners
      aid(i,k)=gsav(j,i)*gsav(j,k)+aid(i,k)
  600 continue
  610 continue
  620 continue
      call sgefa(aid,naid,nindep1,aid2(1,1),info)
      call sgedi(aid,naid,nindep1,aid2(1,1),aid2(1,5),aid2(1,6),0)
      stem=phi
c     signal for srch to read in data
      nnn=ners-nindep1
      if(nnn.ne.0) stem=stem/nnn
      do 660 i=1,nindep1
      do 650 k=1,nindep1
      aid(i,k)=sign(sqrt(stem*abs(aid(i,k))),aid(i,k))
  650 continue
  660 continue
      write(9,64)
      do 670 i=1,nindep1
      write(9,65)  (aid(i,k),k=1,nindep1)
  670 continue
      write(9,59)
      write(9,7)(v(k),k=1,nu)
      gradsq=0.d0
      do 25 k=1,nindep1
   25 gradsq=gradsq+grad(k)*grad(k)
      write(9,9)gradsq
      write(9,10)(grad(k),k=1,nindep1)
      call second(tend)
      tend=tend-tstart
      write(9,54)nprob,tend
      return
  240 if(phi-phimin)250,250,260
c     phi l th or = phimin, conv
  250 write(9,4)
      write(9,13)phi,phimin
      go to 230
  260 phirel=(phiprev-phi)/phiprev
      if(phirel)262,280,265
c     phi gr or = phiprev
  262 write(9,16)phi,phiprev
      go to 230
  265 if(phirel-relphi)270,270,280
c     phirel l th or = relphi, conv
  270 write(9,4)
      write(9,14)phirel,relphi
      go to 230
c     phirel gr th relphi
  280 phiprev=phi
      if(iter-iters)300,300,290
c     iter = iters+1, count exit
  290 write(9,4)
      write(9,15)iters
      if(iter-99999) 230,231,231
c     iter l th or = iters
  300 if(zsq-zsqmin)350,350,360
c     zsq l th or = zsqmin, exit
  350 write(9,4)
      write(9,19)zsq,zsqmin
      go to 230
c     zsq gr  th zsqmin
  360 go to 430
c 360 if(yz)390,380,380
c     yz 0.0 or pos
  380 write(9,21)yz
      go to 230
c     yz neg
  390 if(yz+grdzmin)430,400,400
c     yz convergence
  400 write(9,62)yz,grdzmin
      go to 230
  430 call srch
      write(l,100) (v(i),i=1,nu)
  100 format(1x,*    v(i)=*,
     1 (1x,1pe20.12,e20.12,e20.12,e20.12,e20.12))
      write(l,105) iter
      if(nindep1) 435,235,195
  435 nindep1=-nindep1
      go to 235
    1 format(i4)
    2 format(d8.2)
    4 format('****    ****    ****    ****    ****    ****    ******')
    5 format('iters  =',i4,1x,'number of allowed iterations')
    7 format(1p5e24.14)
    8 format(11h 230, phi =,1pe24.14,22h, parameter vector u =)
    9 format(28h 230, gradient squared mag =,1pe10.2,15h,gradient vec =)
   10 format(1p12e10.2)
   11 format(9h phimin =,1pe10.2,47h, convergence when phi less than non
     czero phimin)
   12 format(9h relphi =,1pe10.2,63h, convergence when relative drop in
     c phi less than input relphi)
   13 format(34h 250, small phi convergence, phi =,1pe24.14,29h is less
     c than input phimin =,1pe24.14)
   14 format(48h 270, convergence, small relative phi decrease =,1pe10.2
     c,28h is less than input relphi =,1pe10.2)
   15 format(54h 290, iteration count exit, iterations = input iters =,i
     c4)
   16 format(* optimiz262, present phi =*,1pe22.14,* is gr th or = previ
     cous phi =*,1pe22.14,/*  accept as convergence only if gradients ar
     ce small*)
   18 format(9h zsqmin =,1pe10.2,53h, exit when squared mag of z vector
     c less than zsqmin)
   19 format(41h small z vector exit, squared mag z vec =,1pe10.2,28h is
     c less than input zsqmin =,1pe10.2)
   20 format(9h sslope =,1pe10.2,93h, exit when negative s-slope of phi-
     ceta curve is more positive than sslope, a negative number)
   21 format(* optimiz 380, error halt, grad.z is positive and =*,1pe10.
     c2)
   31 format(* 370, z vec =*)
   43 format(* opt 410, rel phi drop =*,1pe10.2,* is less than input rel
     cphit =*,1pe10.2,*, switch from w-srch to eta-srch*)
   46 format(* opt 470, error halt, ncstr =*,i3)
   51 format(*1problem*,i3,* begins here*)
   52 format(* opt 180, nuflag(k) nbd(k)    bl(k)    parameter vector v(
     ck)       bu(k)        k*)
   53 format(*               *,i2,*       *,i2,*  *,1pe13.6,* *,1pe22.14
     c,*  *,1pe13.6,*  *,i2)
   54 format(8h problem,i2,16h, running time =,1pe14.6)
   55 format(* parameter flag vec nuflag(k), 0 (1) to vary (hold fixed)*
     c* parameter v(k)*)
   56 format(* parameter bound flag vec nbd(k), 0 (v(k) unbounded), 1 (v
     c(k) has lower bound only), 2 (v(k) bounded below and above)*,/* 3*
     c* (v(k) has upper bound only)*)
   57 format(* parameter lower bound vec bl(k), format(6e12.6),  blank c
     component when v(k) is unbounded below or held fixed*)
   58 format(* parameter upper bound vec bu(k), format(6e12.6),  blank c
     component when v(k) is unbounded above or held fixed*)
   59 format(* opt 230, parameter vec v =*)
   61 format(* grdzmin=*,1pe10.2,*, convergence when grad.z between 0.0*
     c* and -grdzmin*)
   62 format(* optimiz 400, convergence, grad.z=*,1pe10.2,* is abs l th*
     c* grdzmin =*,1pe10.2)
   63 format(* optimiz 420, convergence, previous grad.z =*,1pe13.5,*, p
     cresent grad.z =*,1pe13.5,*, rel change =*,1pe10.2,/* is abs l th r
     celgrdz =*,1pe10.2)
   64 format(1x,*statistical output, standard deviations*)
   65 format(1x,12(x,1pe10.2))
  105 format(*      , end of iteration*,i3,/)
      end




