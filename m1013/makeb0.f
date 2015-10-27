*deck @(#)makeb0.f	5.1  11/6/94
      subroutine makeb0(values,nnp,num,ntriang,isq,ds,b0,c,dlag,
     #                 pt,nindep,alpha,beta,nshell,orbshl,
     #                 minshl,maxshl,numshl,lag,glag,t1,t2,nder,ops)
c
c***begin prologue     makeb0
c***date written       861211  (yymmdd)
c***revision date      930711  (yymmdd)
c   july 7, 1993       rlm at lanl
c      modifying to write the b0 vectors to the tints file as opposed to ints.
c***keywords           b0 vectors, coupled perturbed hartree-fock, cphf
c***author             saxe, paul (lanl)
c***source             @(#)makeb0.f	5.1   11/6/94
c***purpose            form the right hand sides (b0's) of the cphf equations
c***description
c
c     b0(ij,a) = -dlag(ij,a) + dlag(ji,a) +
c
c                all  occ
c                sum  sum ds(kl,a) * { 2 * (alpha(i,k)-alpha(jk)) * [ij;kl] +
c                 k  > l
c
c                         (beta(i,k)-beta(j,k)) * ([ik;jl] + [il;jk]) +
c
c                         delta(kj) * (lag(il)-glag(il,jshell)) -
c
c                         delta(ki) * (lag(jl)-glag(jl,ishell)) }
c
c                                        +
c
c                          occ
c                    0.5 * sum  s(kk,a) * { 2 * (alpha(i,k)-alpha(j,k)) *
c                           k
c
c                               [ij;kk] + 2 * (beta(i,k)-beta(j,k)) * [ik;jk] }
c
c***references         "unified theoretical treatment of analytic first
c                       and second derivatives in open-shell hartree-
c                       fock theory", y. osamura, y. yamaguchi, p. saxe,
c                       m. a. vincent, j. f. gaw and h. f. schaefer iii,
c                       chemical physics 72 (1983) 131-139.
c
c***routines called    iosys
c
c***end prologue       makeb0
c
      implicit integer (a-z)
c
      character*(*) ops
      real*8 values(nnp,ntriang)
      real*8 isq(num,num)
      real*8 alpha(nshell,nshell),beta(nshell,nshell)
      real*8 lag(num,num),glag(num,num,nshell)
      real*8 ds(nnp,nder),b0(nindep,nder),c(num,num)
      real*8 t1(num,num),t2(num,num),dlag(num,num,nder)
      real*8 fac
      integer pt(nnp),orbshl(num)
      integer minshl(nshell),maxshl(nshell),numshl(nshell)
      logical logkey
c
      common /io/ inp,iout
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
      call iosys('read real "scf vector" from rwf',num**2,c,0,' ')
      call iosys('read real "ao derivative overlap integrals" '//
     $     'from rdints',nnp*nder,ds,0,' ')
c
      do 60 der=1,nder
         call trtosq(t1,ds(1,der),num,nnp)
         call ebc(t2,t1,c,num,num,num)
         call ebtc(t1,c,t2,num,num,num)
         call sqtotr(ds(1,der),t1,num,nnp)
   60 continue
c
      call iosys('write real "mo derivative overlap integrals" '//
     $     'to rwf',nnp*nder,ds,0,' ')
c
c     ----- read in the derivative lagrangians and form
c
c             b0(ij,a) = -dlag(i,j,a) + dlag(j,i,a)
c
      call iosys('read real "mo derivative hf lagrangian" from rwf',
     #            num**2*nder,dlag,0,' ')
c
      do 65 ishell=1,nshell
         do 64 jshell=1,ishell-1
            do 63 i=minshl(ishell),maxshl(ishell)
               ia=i*(i-1)/2
               do 62 j=minshl(jshell),maxshl(jshell)
                  pij=pt(ia+j)
                  do 61 der=1,nder
                     b0(pij,der)=-dlag(i,j,der)+dlag(j,i,der)
   61             continue
   62          continue
   63       continue
   64    continue
   65 continue
c
c     ----- rewind the mo integrals -----
c
      call iosys('rewind "mo two-electron integrals" on tints',
     $            0,0,0,' ')
c
c     ----- read through integrals -----
c
      n=0
c
      maxkl=0
      kl=0
      do 110 k=1,num
         k1=k*(k-1)/2+1
         do 100 l=1,k
            l1=l*(l-1)/2+1
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
c           ----- work on coulomb part of b0 -----
c
c             s(kl,a) * 2* (alpha(i,k)-alpha(j,k)) * [ij;kl]
c
            kshell=orbshl(k)
            lshell=orbshl(l)
c
c           ----- orbital k>=l and l must be occupied ------
c
            if (lshell.le.nshell-1) then
               do 4 ishell=2,nshell
                  do 3 jshell=1,ishell-1
                     fac=2.0*(alpha(ishell,kshell)-alpha(jshell,kshell))
                     if (k.eq.l) fac=fac/2.0
c
                     do 2 i=minshl(ishell),maxshl(ishell)
                        ia=i*(i-1)/2
                        do 1 j=minshl(jshell),maxshl(jshell)
                           ij=ia+j
                           pij=pt(ij)
                           do 51 der=1,nder
                              b0(pij,der)=b0(pij,der)+ds(kl,der)*fac*
     #                                    values(ij,kl-minkl+1)
   51                      continue
    1                   continue
    2                continue
    3             continue
    4          continue
            end if
c
c           ----- and exchange part of b0:
c
c            s(kl,a) * (beta(i,k) - beta(j,k)) * ([ik;jl] + [il;jk])
c
c           ----- square up the integrals -----
c
            call trtosq(isq,values(1,kl-minkl+1),num,nnp)
c
c----------------------------------------------------------------------
c  interchanging the indices on the first term (since i label the
c  integral [ij;kl]) gives the following expression:
c
c      s(jl,a) * (beta(i,j)-beta(k,j)) * [ij;kl]
c
c             goes to b0(ik,a)   where
c                    i>k  and j>l
c----------------------------------------------------------------------
c
            do 8 ishell=kshell+1,nshell
               do 7 jshell=lshell,nshell
c
                  fac=beta(ishell,jshell)-beta(kshell,jshell)
c
                  do 6 i=minshl(ishell),maxshl(ishell)
                     ik=i*(i-1)/2+k
                     pik=pt(ik)
c
c                    ----- ensure j>=l and also pick up [ij;kj] terms here -----
c
                     if (jshell.eq.lshell) then
                        minj=l
                     else
                        minj=minshl(jshell)
                     end if
c
                     do 5 j=minj,maxshl(jshell)
                        jl=j*(j-1)/2+l
                        do 52 der=1,nder
                           b0(pik,der)=b0(pik,der)+ds(jl,der)*fac*
     #                                 isq(i,j)
   52                   continue
    5                continue
    6             continue
    7          continue
    8       continue
c
c----------------------------------------------------------------------
c   we also need to flip the k and l indices if k>l since we'll
c   only see [ij;kl] and not [ij;lk]. this gives the following:
c
c      s(jk,a) * (beta(i,j)-beta(l,j)) * [ij;kl]
c
c             goes to b0(il)    where
c                    i>l  and j>=k
c----------------------------------------------------------------------
c
            if (k.ne.l) then
               do 12 ishell=lshell+1,nshell
                  do 11 jshell=kshell,nshell
c
                     fac=beta(ishell,jshell)-beta(lshell,jshell)
c
                     do 10 i=minshl(ishell),maxshl(ishell)
                        il=i*(i-1)/2+l
                        pil=pt(il)
c
c                       ----- ensure j>=k and pick up [ij;lj] terms here
c
                        if (jshell.eq.kshell) then
                           minj=k
                        else
                           minj=minshl(jshell)
                        end if
c
                        do 9 j=minj,maxshl(jshell)
                           jk=j*(j-1)/2+k
                           do 53 der=1,nder
                              b0(pil,der)=b0(pil,der)+ds(jk,der)*fac*
     #                                    isq(i,j)
   53                      continue
    9                   continue
   10                continue
   11             continue
   12          continue
            end if
c
c----------------------------------------------------------------------
c  interchanging the indices on the second term (since i label the
c  integral [ij;kl]) gives the following expression:
c
c       s(lj,a) * (beta(i,l)-beta(k,l)) * [ij;kl]
c
c             goes to b0(ik)    where
c                    i>k  and j<l
c  nb. the j=l case has been done above.
c----------------------------------------------------------------------
c
            do 16 ishell=kshell+1,nshell
               do 15 jshell=1,lshell
c
                  fac=beta(ishell,lshell)-beta(kshell,lshell)
c
                  do 14 i=minshl(ishell),maxshl(ishell)
                     ik=i*(i-1)/2+k
                     pik=pt(ik)
c
c                    ----- ensure j<l -----
c
                     if (jshell.eq.lshell) then
                        jmax=l-1
                     else
                        jmax=maxshl(jshell)
                     end if
c
                     do 13 j=minshl(jshell),jmax
                        lj=l*(l-1)/2+j
                        do 54 der=1,nder
                           b0(pik,der)=b0(pik,der)+ds(lj,der)*fac*
     #                                 isq(i,j)
   54                   continue
   13                continue
   14             continue
   15          continue
   16       continue
c
c----------------------------------------------------------------------
c   we also need to flip the k and l indices if k>l since we'll
c   only see [ij;kl] and not [ij;lk]. this gives the following:
c
c       s(kj,a) * (beta(i,k)-beta(l,k)) * [ij;kl]
c
c             goes to b0(il)    where
c                    i>l  and j<k
c   nb. again, the diagonal term [ij;jl] already accounted for
c----------------------------------------------------------------------
c
            if (k.ne.l) then
               do 20 ishell=lshell+1,nshell
                  do 19 jshell=1,kshell
c
                     fac=beta(ishell,kshell)-beta(lshell,kshell)
c
                     do 18 i=minshl(ishell),maxshl(ishell)
                        il=i*(i-1)/2+l
                        pil=pt(il)
c
c                       ----- ensure j<k -----
c
                        if (jshell.eq.kshell) then
                           jmax=k-1
                        else
                           jmax=maxshl(jshell)
                        end if
c
                        do 17 j=minshl(jshell),jmax
                           kj=k*(k-1)/2+j
                           do 55 der=1,nder
                              b0(pil,der)=b0(pil,der)+ds(kj,der)*fac*
     #                                    isq(i,j)
   55                      continue
   17                   continue
   18                continue
   19             continue
   20          continue
            end if
c
c
  100    continue
  110 continue
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
c            b0(ij) =+ ds(jk,a) * (lag(ik)-glag(i,k,jshell))
c
      do 120 ishell=2,nshell
         do 119 jshell=1,ishell-1
            do 118 kshell=1,jshell
               do 117 i=minshl(ishell),maxshl(ishell)
                  ia=i*(i-1)/2
                  do 116 j=minshl(jshell),maxshl(jshell)
                     ja=j*(j-1)/2
                     pij=pt(ia+j)
c
c                    ----- ensure k<j -----
c
                     if (kshell.eq.jshell) then
                        kmax=j-1
                     else
                        kmax=maxshl(kshell)
                     end if
c
                     do 115 k=minshl(kshell),kmax
                        jk=ja+k
                        do 114 der=1,nder
                           b0(pij,der)=b0(pij,der)+ds(jk,der)*
     #                                   (lag(i,k)-glag(i,k,jshell))
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
c            b0(ij) =- ds(ik,a) * (lag(j,k)-glag(j,k,ishell))
c
      do 130 ishell=2,nshell
         do 129 jshell=1,ishell-1
            do 128 kshell=1,min(ishell,nshell-1)
               do 127 i=minshl(ishell),maxshl(ishell)
                  ia=i*(i-1)/2
                  do 126 j=minshl(jshell),maxshl(jshell)
                     ja=j*(j-1)/2
                     pij=pt(ia+j)
c
c                    ----- ensure k<i -----
c
                     if (kshell.eq.ishell) then
                        kmax=i-1
                     else
                        kmax=maxshl(kshell)
                     end if
c
                     do 125 k=minshl(kshell),kmax
                        ik=ia+k
                        do 124 der=1,nder
                           b0(pij,der)=b0(pij,der)-ds(ik,der)*
     #                                 (lag(j,k)-glag(j,k,ishell))
  124                   continue
  125                continue
  126             continue
  127          continue
  128       continue
  129    continue
  130 continue
c
c     ----- and put the b0 vectors out -----
c
      call iosys('write real b0 to tints',nindep*nder,b0,0,' ')
c
c
      if (logkey(ops,'print=cphf=b0',.false.,' ')) then
         write (iout,140)
 140      format (/,10x,'the b0 vectors:')
          call matout(b0,nindep,nder,nindep,nder,iout)
       end if
c
c
      return
      end
