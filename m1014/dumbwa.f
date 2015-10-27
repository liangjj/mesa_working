*deck @(#)dumbwa.f	5.1  11/6/94
      subroutine dumbwa(values,nnp,num,ntriang,isq,ds,wa,dlag,
     #                 pt,nindep,alpha,beta,nshell,orbshl,
     #                 minshl,maxshl,numshl,lag,glag,t1,t2,nder,
     $     ints,ops)
c
c***begin prologue     dumbwa
c***date written       871114  (yymmdd)
c***revision date      yymmdd  (yymmdd)
c***keywords           wa vectors, second derivatives
c***author             saxe, paul (lanl)
c***source             @(#)dumbwa.f	5.1   11/6/94
c***purpose            form the omega a's (wa's) of the hf second
c                      derivative equations. debug version with all
c                      of the integral held in core
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
c***end prologue       dumbwa
c
      implicit integer (a-z)
c
      character*(*) ops
      logical logkey
      real*8 ints(num,num,num,num)
      real*8 values(nnp,ntriang)
      real*8 isq(num,num)
      real*8 alpha(nshell,nshell),beta(nshell,nshell)
      real*8 lag(num,num),glag(num,num,nshell)
      real*8 ds(nnp,nder),wa(num,num,nder)
      real*8 t1(num,num),t2(num,num),dlag(num,num,nder)
      integer pt(nnp),orbshl(num)
      integer minshl(nshell),maxshl(nshell),numshl(nshell)
c
      common /io/ inp,iout
c
      iadd(i,j)=i*(i-1)/2+j
      madd(i,j)=iadd(max(i,j),min(i,j))
c
c     ----- rewind the mo integrals -----
c
      if (ntriang.ne.nnp) call lnkerr('ntriang.ne.nnp')
      call iosys('read real "mo two-electron integrals" '//
     $     'from tints',nnp**2,values,0,' ')
c
c     ----- square them up
c
      do 304 i=1,num
         do 303 j=1,i
            ij=iadd(i,j)
            call trtosq(ints(1,1,i,j),values(1,ij),num,nnp)
            call trtosq(ints(1,1,j,i),values(1,ij),num,nnp)
 303     continue
 304  continue
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
      call iosys('read real "mo derivative overlap integrals" '//
     $     'from rwf',nnp*nder,ds,0,' ')
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
               wa(i,j,der)=4*dlag(j,i,der)
 61         continue
 64      continue
 65   continue
c
c     ----- two-electron terms -----
c
      do 80 i=1,num
         ishell=orbshl(i)
         do 79 j=1,num
            jshell=orbshl(j)
            do 78 k=1,maxshl(nshell-1)
               kshell=orbshl(k)
               do 77 l=1,maxshl(nshell-1)
                  lshell=orbshl(l)
                  kl=madd(k,l)
                  do 76 der=1,nder
                     wa(i,j,der)=wa(i,j,der)-2*ds(kl,der)*
     $                    (2*alpha(kshell,jshell)*ints(i,j,k,l)+
     $                    beta(kshell,jshell)*(ints(i,k,j,l)+
     $                    ints(i,l,j,k)))
 76               continue
 77            continue
 78         continue
 79      continue
 80   continue
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
c            wa(ij) =- 2 * ds(jk,a) * (lag(ik)-glag(i,k,jshell))
c
      do 120 ishell=1,nshell
         do 119 jshell=1,nshell
            do 118 kshell=1,nshell-1
               do 117 i=minshl(ishell),maxshl(ishell)
                  do 116 j=minshl(jshell),maxshl(jshell)
                     do 115 k=minshl(kshell),maxshl(kshell)
                        jk=madd(j,k)
                        do 114 der=1,nder
                           wa(i,j,der)=wa(i,j,der)-2.0d+00*ds(jk,der)*
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
                           wa(i,j,der)=wa(i,j,der)-4.0d+00*ds(ik,der)*
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
     $     num*num*nder,wa,0,' ')
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
