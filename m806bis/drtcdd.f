*deck @(#)drtcdd.f	1.2  7/30/91
      block data drtcdd
c
      implicit integer (a-z)
c
      character*1 multrf,valenc
      character*3 codes,words*18
c
      common /drtcod/ ncodes,dela(9),delb(9),delele(9)
     #,               ntypes,virtul,occupd,valocc,rescor,resvir,frozen
     #,               valvir,opensh,multi,speshl
      common /drtchr/ codes(9),words(9),multrf,valenc
      common /code/   fzc, fzv, cor, vir, doc, uoc, alp, bet, spe
      common /cases/  casev(9)
c
      data ncodes /9/
      data codes /'fzc','fzv','cor','vir','doc','uoc','alp','bet','spe'/
      data dela  /  0  ,  0  ,  1  ,  0  ,  1  ,  0  ,  0  ,  1  ,  0  /
      data delb  /  0  ,  0  ,  0  ,  0  ,  0  ,  0  ,  1  , -1  ,  0  /
      data delele/  0  ,  0  ,  2  ,  0  ,  2  ,  0  ,  1  ,  1  ,  0  /
      data casev /99999,99999,  4  ,  1  ,  4  ,  1  ,  2  ,  3  ,9999 /
      data ntypes /50/
      data virtul /1/, occupd/4/, valocc/6/, rescor/3/, resvir/2/
     #,    frozen/10/, opensh/8/, multi /7/, valvir/5/, speshl/9/
      data multrf /'/'/
      data valenc /'%'/
      data words  /'frozen core       '
     #,            'frozen virtual    '
     #,            'restricted core   '
     #,            'restricted virtual'
     #,            'doubly occupied   '
     #,            'unoccupied        '
     #,            'alpha occupied    '
     #,            'beta occupied     '
     #,            'special orbital   '/
      data fzc/1/, fzv/2/, cor/3/, vir/4/, doc/5/, uoc/6/, alp/7/
     #,    bet/8/, spe/9/
      end
