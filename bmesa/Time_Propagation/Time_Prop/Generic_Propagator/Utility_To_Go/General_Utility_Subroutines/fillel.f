*deck @(#)fillel.f	5.1  11/6/94
      subroutine fillel(ist,iend,el)
c***begin prologue     fillel
c***date written       850601  (yymmdd)
c***keywords           elements, names
c***revision date      850601  (yymmdd)
c***author             martin, richard (lanl).
c***source
c***purpose            fills an array with the names of the elements.
c***description
c                      module to load the array el with the names of the
c                      elements from atomic number ist to iend.
c                      ist can be zero, in which el starts with 'bq',
c                      or -1, in which case el starts with 'x '.
c
c                      call fillel(ist,iend,el)
c                        ist      starting atomic number.
c                        iend     ending atomic number.
c                        el       character*2 array holding the element names.
c
c***references
c***routines called    (none)
c***end prologue       fillel
      implicit integer(a-z)
      character*(*) el(*)
      character*2 eldat(106),blank
      data maxend/104/,  eldat/'x ','zq','h ','he','li','be','b ','c ',
     $'n ','o ','f ','ne','na','mg','al','si','p ','s ','cl','ar','k ',
     $'ca','sc','ti','v ','cr','mn','fe','co','ni','cu','zn','ga','ge',
     $'as','se','br','kr','rb','sr','y ','zr','nb','mo','tc','ru','rh',
     $'pd','ag','cd','in','sn','sb','te','i ','xe','cs','ba','la','ce',
     $'pr','nd','pm','sm','eu','gd','tb','dy','ho','er','tm','yb','lu',
     $'hf','ta','w ','re','os','ir','pt','au','hg','tl','pb','bi','po',
     $'at','rn','fr','ra','ac','th','pa','u ','np','pu','am','cm','bk',
     $'cf','es','fm','md','no','lr','ky'/, blank/'  '/
      save maxend,eldat,blank
c
c
      ist1 = max(ist,-1)
      iend1 = min(iend,maxend)
      num = iend1 - ist1 + 1
c
      if(num.lt.1) return
      do 10 i = 1, num
   10     el(i) = eldat(i+ist1+1)
c
      num1 = iend - iend1
      if(num1.lt.1) return
      do 20 i = 1, num1
   20     el(i+num) = blank
c
c
      return
      end
