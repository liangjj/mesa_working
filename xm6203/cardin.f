*deck @(#)cardin.f	1.2  10/27/94
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
