*deck @(#)makewa.f	5.1  11/6/94
      subroutine makewa(values,nnp,num,ntriang,isq,ds,wa,dlag,
     #                 pt,nindep,alpha,beta,nshell,orbshl,
     #                 minshl,maxshl,numshl,lag,glag,t1,t2,nder,ops)
c
c***begin prologue     makewa
c***creation date      920805  (yymmdd)
c***revision date      yymmdd  (yymmdd)
c***keywords           wa vectors, second derivatives
c***author             lengsfield, byron (llnl)
c***source             @(#)makewa.f	5.1   11/6/94
c***purpose            form the omega a's (wa's) of the hf second
c                      derivative equations.
c***description
c
c     wa(ij,a) = 4 * dlag(ji,a) -
c
c                occ  occ
c            2 * sum  sum ds(kl,a) * { 2 * alpha(k,j)) * [ij;kl]
c                 k    l
c
c                              + beta(k,j) * ([ik;jl] + [il;jk])}
c
c                         occ
c                  -  2 * sum ds(jk,a) * (lag(ik)+glag(ik,jshell))
c                          k
c
c                         occ
c                  -  4 * sum  s(ik,a) * lag(jk)
c                          k
c
c***references         "unified theoretical treatment of analytic first
c                       and second derivatives in open-shell hartree-
c                       fock theory", y. osamura, y. yamaguchi, p. saxe,
c                       m. a. vincent, j. f. gaw and h. f. schaefer iii,
c                       chemical physics 72 (1983) 131-139.
c
c***routines called    iosys
c
c***end prologue       makewa
c
      implicit integer (a-z)
c
      character*(*) ops
      logical logkey
      real*8 values(nnp,ntriang)
      real*8 isq(num,num)
      real*8 alpha(nshell,nshell),beta(nshell,nshell)
      real*8 lag(num,num),glag(num,num,nshell)
      real*8 ds(nnp,nder),wa(num,num,nder)
      real*8 t1(num,num),t2(num,num),dlag(num,num,nder)
      real*8 fac,ff,fx,two,four
      integer pt(nnp),orbshl(num)
      integer minshl(nshell),maxshl(nshell),numshl(nshell)
c
      parameter (two=2.0d+00,four=4.0d+00)
c
      common /io/ inp,iout
c
      iadd(i,j)=i*(i-1)/2+j
      madd(i,j)=iadd(max(i,j),min(i,j))
c
c     ----- create the array of orbital shells (shells) -----
c
      do 58 ishell=1,nshell
         do 57 i=minshl(ishell),maxshl(ishell)
            orbshl(i)=ishell
   57    continue
   58 continue
c
c     ----- read in and transform the derivative overlap integrals -----
c
c..bhl
      call iosys('read real "mo derivative overlap integrals" '//
     $           'from rwf',nnp*nder,ds,0,' ')
c..bhl
c
c     ----- read in the derivative lagrangians and form
c
c             wa(ij,a) = 4 * dlag(j,i,a)
c
      call iosys('read real "mo derivative hf lagrangian" from rwf',
     #            num**2*nder,dlag,0,' ')
c
      do 65 i=1,num
         do 64 j=1,num
            do 61 der=1,nder
               wa(i,j,der)=four*dlag(j,i,der)
 61         continue
 64      continue
 65   continue
c
c     ----- rewind the mo integrals -----
c
      call iosys('rewind "mo two-electron integrals" on tints',
     $           0,0,0,' ')
c
c     ----- read through integrals -----
c
      n=0
c
      noc=maxshl(nshell-1)
c
c     write(iout,*)' noc ',noc
c
      maxkl=0
      kl=0
      do 110 k=1,noc
         do 100 l=1,k
            kl=kl+1
c
c     ----- check that this triangle of integrals is in core -----
c
            if (kl.gt.maxkl) then
               minkl=maxkl+1
               maxkl=min(nnp,maxkl+ntriang)
               lnread=(maxkl-minkl+1)*nnp
               call iosys('read real "mo two-electron integrals" '//
     $              'from tints without rewinding',lnread,values,0,' ')
            end if
c
c           ----- square up the integrals -----
c
            call trtosq(isq,values(1,kl-minkl+1),num,nnp)
c
c           ----- work on coulomb part of wa -----
c
c             s(kl,a) * 4 * alpha(k,j) * [ij;kl]
c
            kshell=orbshl(k)
            lshell=orbshl(l)
c
c           ----- orbital k>=l and k,l must be occupied ------
c               take care of kl and lk simultaneously, hence
c               two 'alpha' terms.
c
            if(k.eq.l) then
               fac=-two
            else
               fac=-four
            end if
c
            do 4 jshell=1,nshell
               ff=fac*(alpha(kshell,jshell)+alpha(lshell,jshell))
               js=minshl(jshell)
               jn=maxshl(jshell)-js+1
               len=jn*num
               do 3 der=1,nder
                  fx=ff*ds(kl,der)
                  call saxpy(len,fx,isq(1,js),1,wa(1,js,der),1) 
    3          continue
    4       continue
c
c           ----- and exchange part of wa:
c
c            -2 * s(kl,a) * beta(k,j) * ([ik;jl] + [il;jk])
c
c----------------------------------------------------------------------
c
c      -2 * s(jl,a) * beta(j,k) * [ij;kl]
c
c----------------------------------------------------------------------
c
            do 8 jshell=1,nshell
               fac=-two*beta(lshell,jshell)
               do 7 j=minshl(jshell),maxshl(jshell)
                  do 6 i=1,noc
                     t1(j,i)=fac*isq(j,i)
  6               continue
  7            continue
  8         continue
c
            do 52 der=1,nder
               do 51 i=1,noc
                  li=madd(l,i)
                  call saxpy(num,ds(li,der),t1(1,i),1,wa(k,1,der),num)   
 51            continue
 52         continue
c
            if(k.ne.l) then
c
c   flip    k <-> l  indices
c
               do 88 jshell=1,nshell
                  fac=-two*beta(kshell,jshell)
                  do 77 j=minshl(jshell),maxshl(jshell)
                     do 66 i=1,noc
                        t1(j,i)=fac*isq(j,i)
 66                  continue
 77               continue
 88            continue
c
               do 522 der=1,nder
                  do 511 i=1,noc
                     ki=madd(k,i)
                     call saxpy(num,ds(ki,der),t1(1,i),1,
     $                          wa(l,1,der),num)   
 511              continue
 522           continue
            end if
c
c----------------------------------------------------------------------
c
c       -2 * s(lj,a) * beta(l,k) * [ij;kl]
c
c----------------------------------------------------------------------
c
            fac=-two*beta(lshell,kshell)
            do 11 j=1,noc
               do 12 i=1,num
                  t1(i,j)=fac*isq(i,j)
 12            continue
 11         continue
c
            do 13 der=1,nder
               do 14 j=1,noc
                  lj=madd(l,j)
                  call saxpy(num,ds(lj,der),t1(1,j),1,wa(1,k,der),1)
 14            continue
 13         continue
c
c
c           ----- flip   k <-> l -----
            if(k.ne.l) then
               do 133 der=1,nder
                  do 144 j=1,noc
                     kj=madd(k,j)
                     call saxpy(num,ds(kj,der),t1(1,j),1,wa(1,l,der),1)
 144              continue
 133           continue
            end if
c
  100    continue
  110 continue
c
c     ----- finish exchange terms -----
c
      do 210 k=noc+1,num
         kshell=orbshl(k)
         do 200 l=1,noc
            kl=kl+1
c
            lshell=orbshl(l)
c
c     ----- check that this triangle of integrals is in core -----
c
            if (kl.gt.maxkl) then
               minkl=maxkl+1
               maxkl=min(nnp,maxkl+ntriang)
               lnread=(maxkl-minkl+1)*nnp
               call iosys('read real "mo two-electron integrals" '//
     $              'from tints without rewinding',lnread,values,0,' ')
            end if
c
c           ----- square up the integrals -----
c
            call trtosq(isq,values(1,kl-minkl+1),num,nnp)
c
c----------------------------------------------------------------------
c
c      -2 * s(jl,a) * beta(j,k) * [ij;kl]
c
c----------------------------------------------------------------------
c
            do 89 jshell=1,nshell
               fac=-two*beta(lshell,jshell)
               do 79 j=minshl(jshell),maxshl(jshell)
                  do 69 i=1,noc
                     t1(j,i)=fac*isq(j,i)
  69              continue
  79           continue
  89        continue
c
            do 72 der=1,nder
               do 71 i=1,noc
                  li=madd(l,i)
                  call saxpy(num,ds(li,der),t1(1,i),1,wa(k,1,der),num)   
 71            continue
 72         continue
c
c----------------------------------------------------------------------
c
c       -2 * s(lj,a) * beta(l,k) * [ij;kl]
c
c----------------------------------------------------------------------
c
            fac=-two*beta(lshell,kshell)
            do 81 j=1,noc
               do 82 i=1,num
                  t1(i,j)=fac*isq(i,j)
 82            continue
 81         continue
c
            do 83 der=1,nder
               do 84 j=1,noc
                  lj=madd(l,j) 
                  call saxpy(num,ds(lj,der),t1(1,j),1,wa(1,k,der),1)
 84            continue
 83         continue
c
c
 200     continue
c
         do 300 l=noc+1,k
            kl=kl+1
c
c     ----- check that this triangle of integrals is in core -----
c
            if (kl.gt.maxkl) then
               minkl=maxkl+1
               maxkl=min(nnp,maxkl+ntriang)
               lnread=(maxkl-minkl+1)*nnp
               call iosys('read real "mo two-electron integrals" '//
     $              'from tints without rewinding',lnread,values,0,' ')
            end if
 300     continue
c
 210  continue
c
c           ----- square up the integrals -----
c


c
c     ----- get the lagrangians ------
c
      call rzero(glag,num**2*nshell)
      call iosys('read real "scf mo lagrangian" from rwf',num**2,lag,0,
     #           ' ')
      call iosys('read real "generalised hf lagrangian" from rwf',
     #            num**2*(nshell-1),glag,0,' ')
c
c     ----- add in the first lagrangian term:
c
c            wa(ij) =- 2 * ds(jk,a) * (lag(ik)+glag(i,k,jshell))
c
      do 120 ishell=1,nshell
         do 119 jshell=1,nshell
            do 118 kshell=1,nshell-1
               do 117 i=minshl(ishell),maxshl(ishell)
                  do 116 j=minshl(jshell),maxshl(jshell)
                     do 115 k=minshl(kshell),maxshl(kshell)
                        jk=madd(j,k)
                        do 114 der=1,nder
                           wa(i,j,der)=wa(i,j,der)-two*ds(jk,der)*
     #                                   (lag(i,k)+glag(i,k,jshell))
  114                   continue
  115                continue
  116             continue
  117          continue
  118       continue
  119    continue
  120 continue
c
c     ----- add in the second lagrangian term:
c
c            wa(ij) =- 4 * ds(ik,a) * lag(j,k)
c
      do 130 ishell=1,nshell
         do 129 jshell=1,nshell
            do 128 kshell=1,nshell-1
               do 127 i=minshl(ishell),maxshl(ishell)
                  do 126 j=minshl(jshell),maxshl(jshell)
                     do 125 k=minshl(kshell),maxshl(kshell)
                        ik=madd(i,k)
                        do 124 der=1,nder
                           wa(i,j,der)=wa(i,j,der)-four*ds(ik,der)*
     #                                 lag(j,k)
  124                   continue
  125                continue
  126             continue
  127          continue
  128       continue
  129    continue
  130 continue
c
c     ----- and put the wa vectors out -----
c
      call iosys('write real "derivative omegas" to rwf',
     $           num*num*nder,wa,0,' ')
c
c     ----- and print them -----
c
      if (logkey(ops,'print=gradient=omega',.false.,' ')) then
         do 220 der=1,nder
            write (iout,219) der
 219        format (5x,'omega, derivative ',i3)
            call matout(wa(1,1,der),num,num,num,num,iout)
 220     continue
      end if
c
c
      return
      end
