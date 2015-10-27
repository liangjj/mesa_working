*deck @(#)trans2.f	5.1  11/6/94
      subroutine trans2(orig,half,conint,ia,ib,ic,id,fa,fb,fc,fd,
     #                  a,b,c,d,t1,t2,lenblk)
c
c***module to transform two-electron integrals from the primitive-ao
c   basis to contracted functions, vectorized over the number of
c   integrals in the shell-block.
c
c paul saxe                6 august 1984                      lanl
c
      implicit integer (a-z)
c
      real*8 orig(ia*ib*ic*id*lenblk),half(ic*id*lenblk*fa*fb)
      real*8 conint(lenblk*fa*fb*fc*fd),a(ia,fa),b(ib,fb),c(ic,fc)
      real*8 d(id,fd),t1(fa*ib*ic*id*lenblk),t2(fc*id*lenblk*fa*fb)
      real*8 small
c
      data small /1.0d-09/
      save small
c
c     ----- timing -----
c
c
c     ----- zero accumulator for half-transformed integrals -----
c
      call rzero(half,ic*id*fa*fb*lenblk)
c
c     ----- do first half-transformation -----
c
c
c     ----- first quarter-transformation -----
c
      call rzero(t1,fa*ib*ic*id*lenblk)
      do 2 aa=1,fa
         do 1 k=1,ia
            if (abs(a(k,aa)).lt.small) go to 1
            call saxpy(ib*ic*id*lenblk,a(k,aa),orig(k),ia,t1(aa),fa)
c            call gsaxpy(t1(aa),fa,orig(k),ia,a(k,aa),ib*ic*id*lenblk)
    1    continue
    2 continue
c
c     ----- second quarter-transformation -----
c
      do 6 bb=1,fb
         do 5 aa=1,fa
            do 4 k=1,ib
               if (abs(b(k,bb)).lt.small) go to 4
               call saxpy(ic*id*lenblk,b(k,bb),t1(aa+fa*(k-1)),fa*ib,
     #                    half(1+ic*id*lenblk*((aa-1)+fa*(bb-1))),1)
c               call gsaxpy(half(1+ic*id*lenblk*((aa-1)+fa*(bb-1))),1,
c     #                     t1(aa+fa*(k-1)),fa*ib,b(k,bb),ic*id*lenblk)
    4       continue
    5    continue
    6 continue
c
c     ----- second half-transformation -----
c
      call rzero(conint,fa*fb*fc*fd*lenblk)
c
c
c     ----- third quarter-transformation -----
c
      call rzero(t2,fc*id*lenblk*fa*fb)
      do 10 cc=1,fc
         do 9 k=1,ic
            if (abs(c(k,cc)).lt.small) go to 9
            call saxpy(id*lenblk*fa*fb,c(k,cc),half(k),ic,t2(cc),fc)
c            call gsaxpy(t2(cc),fc,half(k),ic,c(k,cc),id*lenblk*fa*fb)
    9    continue
   10 continue
c
c     ----- fourth quarter-transformation -----
c
      do 14 dd=1,fd
         do 13 cc=1,fc
            do 12 k=1,id
               if (abs(d(k,dd)).lt.small) go to 12
               call saxpy(lenblk*fa*fb,d(k,dd),t2(cc+fc*(k-1)),fc*id,
     #                    conint(1+lenblk*fa*fb*((cc-1)+fc*(dd-1))),1)
c               call gsaxpy(conint(1+lenblk*fa*fb*((cc-1)+fc*(dd-1))),1,
c     #                     t2(cc+fc*(k-1)),fc*id,d(k,dd),lenblk*fa*fb)
   12       continue
   13    continue
   14 continue
c
c     ----- timing -----
c
c
c
      return
      end
