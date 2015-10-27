h48787
s 00021/00000/00000
d D 1.1 94/02/16 20:34:56 mesa 1 0
c date and time created 94/02/16 20:34:56 by mesa
e
u
U
f e 0
t
T
I 1
*deck %W%  %G%
c----------------------------------------------------------------------c
c                 card input routine                                   c
c                 reads from first card until a $  encountered         c
c----------------------------------------------------------------------c
      subroutine cardin (card)
      common /io/ inp, iout
      character *(*) card
      ist=1
      ifin=80
   10 read (inp,30) card(ist:ifin)
      if (card(ist:ist).eq.'$') go to 20
      if (card(ist:ist).eq.'%') go to 20
      if (card(ist:ist).eq.'#') go to 20
      ist=ifin+1
      ifin=ifin+80
      go to 10
   20 return
c
   30 format (a80)
      end
E 1
