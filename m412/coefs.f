*deck @(#)coefs.f	5.1  11/6/94
      subroutine coefs(ex,cf,ishell,typeno,typnam,ntypes)
c***begin prologue     coefs
c***date written       850601  yymmdd
c***revision date      yymmdd  yymmdd
c***keywords           sto-3g, basis sets, initial guess
c***author             martin, richard (lanl)
c***source             @(#)coefs.f	5.1   11/6/94
c***purpose            loads sto-3g exponents and contraction coefficients.
c***description
c     call coefs(ex,cf,ishell,typeno,typnam,ntypes)
c       ex      exponent array(3).
c       cf      contraction coefficient array(3).
c       ishell  g82 orbital type:
c                 1  1s
c                 2  2s
c       typeno  orbital type.
c       typnam  mesa shell type('s','p',etc.).
c       ntypes  number of shell types.
c***references
c***routines called    (none)
c***end prologue       coefs
      implicit integer(a-z)
      real*8 ex(3),cf(3)
      character*8 type
      character*(*) typnam(ntypes)
c
c     sto-3g exponents and contraction coefficients.
c     ishell denotes the orbital type.
c     note that 6s, 6p, 7s, 7p, and 6d use lower expansions
c     at the moment.
c
      if(ishell.lt.1.or.ishell.gt.19)
     $   call lnkerr(' unrecognized shell in coefs')
      goto(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19), ishell
c
c     1s
    1 ex(1)=2.227660584d00
      cf(1)=1.543289673d-01
      ex(2)=4.057711562d-01
      cf(2)=5.353281423d-01
      ex(3)=1.098175104d-01
      cf(3)=4.446345422d-01
      type='s'
      goto 100
c
c     2s
    2 ex(1)=9.942027296d-01
      cf(1)=-9.996722919d-02
      ex(2)=2.310313333d-01
      cf(2)=3.995128261d-01
      ex(3)=7.513856000d-02
      cf(3)=7.001154689d-01
      type='s'
      goto 100
c
c     2p
    3 ex(1)=9.942027296d-01
      cf(1)=1.559162750d-01
      ex(2)=2.310313333d-01
      cf(2)=6.076837186d-01
      ex(3)=7.513856000d-02
      cf(3)=3.919573931d-01
      type='p'
      goto 100
c
c     3s
    4 ex(1)=4.828540806d-01
      cf(1)=-2.196203690d-01
      ex(2)=1.347150629d-01
      cf(2)=2.255954336d-01
      ex(3)=5.272656258d-02
      cf(3)=9.003984260d-01
      type='s'
      goto 100
c
c     3p
    5 ex(1)=4.828540806d-01
      cf(1)=1.058760429d-02
      ex(2)=1.347150629d-01
      cf(2)=5.951670053d-01
      ex(3)=5.272656258d-02
      cf(3)=4.620010120d-01
      type='p'
      goto 100
c
c     4s
    6 ex(1)= 2.464581400d-01
      cf(1)=-3.088441215d-01
      ex(2)= 9.095855374d-02
      cf(2)= 1.960641166d-02
      ex(3)= 4.016825636d-02
      cf(3)= 1.131034442d+00
      type='s'
      goto 100
c
c     4p
    7 ex(1)= 2.464581400d-01
      cf(1)=-1.215468600d-01
      ex(2)= 9.095855374d-02
      cf(2)= 5.715227604d-01
      ex(3)= 4.016825636d-02
      cf(3)= 5.498949471d-01
      type='p'
      goto 100
c
c     3d
    8 ex(1)= 5.229112225d-01
      cf(1)= 1.686596060d-01
      ex(2)= 1.639595876d-01
      cf(2)= 5.847984817d-01
      ex(3)= 6.386630021d-02
      cf(3)= 4.056779523d-01
      type='d'
      goto 100
c
c     5s
    9 ex(1)=1.080198458d-01
      cf(1)=-6.617401158d-01
      ex(2)=4.408119382d-02
      cf(2)= 7.467595004d-01
      ex(3)=2.610811810d-02
      cf(3)= 7.146490945d-01
      type='s'
      goto 100
c
c     5p
   10 ex(1)=2.127482317d-01
      cf(1)=-1.389529695d-01
      ex(2)=4.729648620d-02
      cf(2)= 8.076691064d-01
      ex(3)=2.604865324d-02
      cf(3)= 2.726029342d-01
      type='p'
      goto 100
c
c     4d
   11 ex(1)=1.777717219d-01
      cf(1)= 2.308552718d-01
      ex(2)=8.040647350d-02
      cf(2)= 6.042409177d-01
      ex(3)=3.949855551d-02
      cf(3)= 2.595768926d-01
      type='d'
      goto 100
c
c     6s
   12 continue
      goto 9
c
c     6p
   13 continue
      goto 10
c
c     5d
   14 ex(1)=4.913352950d-01
      cf(1)=-2.010175008d-02
      ex(2)=7.329090601d-02
      cf(2)=5.899370608d-01
      ex(3)=3.594209290d-02
      cf(3)=4.658445960d-01
      type='d'
      goto 100
c
c     4f
   15 ex(1)=3.483826963d-01
      cf(1)= 1.737856685d-01
      ex(2)=1.249380537d-01
      cf(2)= 5.973380628d-01
      ex(3)=5.349995725d-02
      cf(3)= 3.929395614d-01
      type='f'
      goto 100
c
c     7s
   16 continue
cpws      goto 9
c   thorium basis
      ex(1)=0.1628
      cf(1)=-1.4073459
      ex(2)=0.09045
      cf(2)=2.2876533
      ex(3)=0.02721
      cf(3)=-0.0964597
      type='s'
      goto 100
c
c     7p
   17 continue
cpws      goto 10
c th basis
      ex(1)=1.135
      cf(1)=-0.3399996
      ex(2)=0.5368
      cf(2)=0.8883881
      ex(3)=0.2106
      cf(3)=0.4111816
      type='p'
      goto 100
c
c     6d
   18 continue
cpws      goto 14
c th basis
      ex(1)=0.3249
      cf(1)=0.6633724
      ex(2)=0.1185
      cf(2)=0.4317079
      ex(3)=0.03
      cf(3)=0.01
      type='d'
      goto 100
c
c     5f
   19 continue
      ex(1)=1.209
      cf(1)=0.5109002
      ex(2)=0.3969
      cf(2)=0.4563037
      ex(3)=0.1205
      cf(3)=0.104856
cpws   19 ex(1)=1.649233885d-01
cpws      cf(1)=1.909729355d-01
cpws      ex(2)=7.487066646d-02
cpws      cf(2)=6.146060459d-01
cpws      ex(3)=3.735787219d-02
cpws      cf(3)=3.059611271d-01
      type='f'
      goto 100
  100 continue
c
c
      do 200 i=1,ntypes
         if(type.eq.typnam(i)) typeno=i
  200 continue
c
c
      return
      end
