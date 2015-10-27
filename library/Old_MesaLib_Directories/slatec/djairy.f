*deck djairy
      subroutine djairy (x, rx, c, ai, dai)
c***begin prologue  djairy
c***subsidiary
c***purpose  subsidiary to dbesj and dbesy
c***library   slatec
c***type      double precision (jairy-s, djairy-d)
c***author  amos, d. e., (snla)
c           daniel, s. l., (snla)
c           weston, m. k., (snla)
c***description
c
c                  djairy computes the airy function ai(x)
c                   and its derivative dai(x) for dasyjy
c
c                                   input
c
c         x - argument, computed by dasyjy, x unrestricted
c        rx - rx=sqrt(abs(x)), computed by dasyjy
c         c - c=2.*(abs(x)**1.5)/3., computed by dasyjy
c
c                                  output
c
c        ai - value of function ai(x)
c       dai - value of the derivative dai(x)
c
c***see also  dbesj, dbesy
c***routines called  (none)
c***revision history  (yymmdd)
c   750101  date written
c   890531  changed all specific intrinsics to generic.  (wrb)
c   891009  removed unreferenced variable.  (wrb)
c   891214  prologue converted to version 4.0 format.  (bab)
c   900328  added type section.  (wrb)
c   910408  updated the author section.  (wrb)
c***end prologue  djairy
c
      integer i, j, m1, m1d, m2, m2d, m3, m3d, m4, m4d, n1, n1d, n2,
     1 n2d, n3, n3d, n4, n4d
      double precision a,ai,ajn,ajp,ak1,ak2,ak3,b,c,ccv,con2,
     1 con3, con4, con5, cv, da, dai, dajn, dajp, dak1, dak2, dak3,
     2 db, ec, e1, e2, fpi12, f1, f2, rtrx, rx, scv, t, temp1, temp2,
     3 tt, x
      dimension ajp(19), ajn(19), a(15), b(15)
      dimension ak1(14), ak2(23), ak3(14)
      dimension dajp(19), dajn(19), da(15), db(15)
      dimension dak1(14), dak2(24), dak3(14)
      save n1, n2, n3, n4, m1, m2, m3, m4, fpi12, con2, con3,
     1 con4, con5, ak1, ak2, ak3, ajp, ajn, a, b,
     2 n1d, n2d, n3d, n4d, m1d, m2d, m3d, m4d, dak1, dak2, dak3,
     3 dajp, dajn, da, db
      data n1,n2,n3,n4/14,23,19,15/
      data m1,m2,m3,m4/12,21,17,13/
      data fpi12,con2,con3,con4,con5/
     1 1.30899693899575d+00, 5.03154716196777d+00, 3.80004589867293d-01,
     2 8.33333333333333d-01, 8.66025403784439d-01/
      data ak1(1), ak1(2), ak1(3), ak1(4), ak1(5), ak1(6), ak1(7),
     1     ak1(8), ak1(9), ak1(10),ak1(11),ak1(12),ak1(13),
     2     ak1(14)         / 2.20423090987793d-01,-1.25290242787700d-01,
     3 1.03881163359194d-02, 8.22844152006343d-04,-2.34614345891226d-04,
     4 1.63824280172116d-05, 3.06902589573189d-07,-1.29621999359332d-07,
     5 8.22908158823668d-09, 1.53963968623298d-11,-3.39165465615682d-11,
     6 2.03253257423626d-12,-1.10679546097884d-14,-5.16169497785080d-15/
      data ak2(1), ak2(2), ak2(3), ak2(4), ak2(5), ak2(6), ak2(7),
     1     ak2(8), ak2(9), ak2(10),ak2(11),ak2(12),ak2(13),ak2(14),
     2     ak2(15),ak2(16),ak2(17),ak2(18),ak2(19),ak2(20),ak2(21),
     3     ak2(22),ak2(23) / 2.74366150869598d-01, 5.39790969736903d-03,
     4-1.57339220621190d-03, 4.27427528248750d-04,-1.12124917399925d-04,
     5 2.88763171318904d-05,-7.36804225370554d-06, 1.87290209741024d-06,
     6-4.75892793962291d-07, 1.21130416955909d-07,-3.09245374270614d-08,
     7 7.92454705282654d-09,-2.03902447167914d-09, 5.26863056595742d-10,
     8-1.36704767639569d-10, 3.56141039013708d-11,-9.31388296548430d-12,
     9 2.44464450473635d-12,-6.43840261990955d-13, 1.70106030559349d-13,
     1-4.50760104503281d-14, 1.19774799164811d-14,-3.19077040865066d-15/
      data ak3(1), ak3(2), ak3(3), ak3(4), ak3(5), ak3(6), ak3(7),
     1     ak3(8), ak3(9), ak3(10),ak3(11),ak3(12),ak3(13),
     2     ak3(14)         / 2.80271447340791d-01,-1.78127042844379d-03,
     3 4.03422579628999d-05,-1.63249965269003d-06, 9.21181482476768d-08,
     4-6.52294330229155d-09, 5.47138404576546d-10,-5.24408251800260d-11,
     5 5.60477904117209d-12,-6.56375244639313d-13, 8.31285761966247d-14,
     6-1.12705134691063d-14, 1.62267976598129d-15,-2.46480324312426d-16/
      data ajp(1), ajp(2), ajp(3), ajp(4), ajp(5), ajp(6), ajp(7),
     1     ajp(8), ajp(9), ajp(10),ajp(11),ajp(12),ajp(13),ajp(14),
     2     ajp(15),ajp(16),ajp(17),ajp(18),
     3     ajp(19)         / 7.78952966437581d-02,-1.84356363456801d-01,
     4 3.01412605216174d-02, 3.05342724277608d-02,-4.95424702513079d-03,
     5-1.72749552563952d-03, 2.43137637839190d-04, 5.04564777517082d-05,
     6-6.16316582695208d-06,-9.03986745510768d-07, 9.70243778355884d-08,
     7 1.09639453305205d-08,-1.04716330588766d-09,-9.60359441344646d-11,
     8 8.25358789454134d-12, 6.36123439018768d-13,-4.96629614116015d-14,
     9-3.29810288929615d-15, 2.35798252031104d-16/
      data ajn(1), ajn(2), ajn(3), ajn(4), ajn(5), ajn(6), ajn(7),
     1     ajn(8), ajn(9), ajn(10),ajn(11),ajn(12),ajn(13),ajn(14),
     2     ajn(15),ajn(16),ajn(17),ajn(18),
     3     ajn(19)         / 3.80497887617242d-02,-2.45319541845546d-01,
     4 1.65820623702696d-01, 7.49330045818789d-02,-2.63476288106641d-02,
     5-5.92535597304981d-03, 1.44744409589804d-03, 2.18311831322215d-04,
     6-4.10662077680304d-05,-4.66874994171766d-06, 7.15218807277160d-07,
     7 6.52964770854633d-08,-8.44284027565946d-09,-6.44186158976978d-10,
     8 7.20802286505285d-11, 4.72465431717846d-12,-4.66022632547045d-13,
     9-2.67762710389189d-14, 2.36161316570019d-15/
      data a(1),   a(2),   a(3),   a(4),   a(5),   a(6),   a(7),
     1     a(8),   a(9),   a(10),  a(11),  a(12),  a(13),  a(14),
     2     a(15)           / 4.90275424742791d-01, 1.57647277946204d-03,
     3-9.66195963140306d-05, 1.35916080268815d-07, 2.98157342654859d-07,
     4-1.86824767559979d-08,-1.03685737667141d-09, 3.28660818434328d-10,
     5-2.57091410632780d-11,-2.32357655300677d-12, 9.57523279048255d-13,
     6-1.20340828049719d-13,-2.90907716770715d-15, 4.55656454580149d-15,
     7-9.99003874810259d-16/
      data b(1),   b(2),   b(3),   b(4),   b(5),   b(6),   b(7),
     1     b(8),   b(9),   b(10),  b(11),  b(12),  b(13),  b(14),
     2     b(15)           / 2.78593552803079d-01,-3.52915691882584d-03,
     3-2.31149677384994d-05, 4.71317842263560d-06,-1.12415907931333d-07,
     4-2.00100301184339d-08, 2.60948075302193d-09,-3.55098136101216d-11,
     5-3.50849978423875d-11, 5.83007187954202d-12,-2.04644828753326d-13,
     6-1.10529179476742d-13, 2.87724778038775d-14,-2.88205111009939d-15,
     7-3.32656311696166d-16/
      data n1d,n2d,n3d,n4d/14,24,19,15/
      data m1d,m2d,m3d,m4d/12,22,17,13/
      data dak1(1), dak1(2), dak1(3), dak1(4), dak1(5), dak1(6),
     1     dak1(7), dak1(8), dak1(9), dak1(10),dak1(11),dak1(12),
     2    dak1(13),dak1(14)/ 2.04567842307887d-01,-6.61322739905664d-02,
     3-8.49845800989287d-03, 3.12183491556289d-03,-2.70016489829432d-04,
     4-6.35636298679387d-06, 3.02397712409509d-06,-2.18311195330088d-07,
     5-5.36194289332826d-10, 1.13098035622310d-09,-7.43023834629073d-11,
     6 4.28804170826891d-13, 2.23810925754539d-13,-1.39140135641182d-14/
      data dak2(1), dak2(2), dak2(3), dak2(4), dak2(5), dak2(6),
     1     dak2(7), dak2(8), dak2(9), dak2(10),dak2(11),dak2(12),
     2     dak2(13),dak2(14),dak2(15),dak2(16),dak2(17),dak2(18),
     3     dak2(19),dak2(20),dak2(21),dak2(22),dak2(23),
     4     dak2(24)        / 2.93332343883230d-01,-8.06196784743112d-03,
     5 2.42540172333140d-03,-6.82297548850235d-04, 1.85786427751181d-04,
     6-4.97457447684059d-05, 1.32090681239497d-05,-3.49528240444943d-06,
     7 9.24362451078835d-07,-2.44732671521867d-07, 6.49307837648910d-08,
     8-1.72717621501538d-08, 4.60725763604656d-09,-1.23249055291550d-09,
     9 3.30620409488102d-10,-8.89252099772401d-11, 2.39773319878298d-11,
     1-6.48013921153450d-12, 1.75510132023731d-12,-4.76303829833637d-13,
     2 1.29498241100810d-13,-3.52679622210430d-14, 9.62005151585923d-15,
     3-2.62786914342292d-15/
      data dak3(1), dak3(2), dak3(3), dak3(4), dak3(5), dak3(6),
     1     dak3(7), dak3(8), dak3(9), dak3(10),dak3(11),dak3(12),
     2    dak3(13),dak3(14)/ 2.84675828811349d-01, 2.53073072619080d-03,
     3-4.83481130337976d-05, 1.84907283946343d-06,-1.01418491178576d-07,
     4 7.05925634457153d-09,-5.85325291400382d-10, 5.56357688831339d-11,
     5-5.90889094779500d-12, 6.88574353784436d-13,-8.68588256452194d-14,
     6 1.17374762617213d-14,-1.68523146510923d-15, 2.55374773097056d-16/
      data dajp(1), dajp(2), dajp(3), dajp(4), dajp(5), dajp(6),
     1     dajp(7), dajp(8), dajp(9), dajp(10),dajp(11),dajp(12),
     2     dajp(13),dajp(14),dajp(15),dajp(16),dajp(17),dajp(18),
     3     dajp(19)        / 6.53219131311457d-02,-1.20262933688823d-01,
     4 9.78010236263823d-03, 1.67948429230505d-02,-1.97146140182132d-03,
     5-8.45560295098867d-04, 9.42889620701976d-05, 2.25827860945475d-05,
     6-2.29067870915987d-06,-3.76343991136919d-07, 3.45663933559565d-08,
     7 4.29611332003007d-09,-3.58673691214989d-10,-3.57245881361895d-11,
     8 2.72696091066336d-12, 2.26120653095771d-13,-1.58763205238303d-14,
     9-1.12604374485125d-15, 7.31327529515367d-17/
      data dajn(1), dajn(2), dajn(3), dajn(4), dajn(5), dajn(6),
     1     dajn(7), dajn(8), dajn(9), dajn(10),dajn(11),dajn(12),
     2     dajn(13),dajn(14),dajn(15),dajn(16),dajn(17),dajn(18),
     3     dajn(19)        / 1.08594539632967d-02, 8.53313194857091d-02,
     4-3.15277068113058d-01,-8.78420725294257d-02, 5.53251906976048d-02,
     5 9.41674060503241d-03,-3.32187026018996d-03,-4.11157343156826d-04,
     6 1.01297326891346d-04, 9.87633682208396d-06,-1.87312969812393d-06,
     7-1.50798500131468d-07, 2.32687669525394d-08, 1.59599917419225d-09,
     8-2.07665922668385d-10,-1.24103350500302d-11, 1.39631765331043d-12,
     9 7.39400971155740d-14,-7.32887475627500d-15/
      data da(1),  da(2),  da(3),  da(4),  da(5),  da(6),  da(7),
     1     da(8),  da(9),  da(10), da(11), da(12), da(13), da(14),
     2     da(15)          / 4.91627321104601d-01, 3.11164930427489d-03,
     3 8.23140762854081d-05,-4.61769776172142d-06,-6.13158880534626d-08,
     4 2.87295804656520d-08,-1.81959715372117d-09,-1.44752826642035d-10,
     5 4.53724043420422d-11,-3.99655065847223d-12,-3.24089119830323d-13,
     6 1.62098952568741d-13,-2.40765247974057d-14, 1.69384811284491d-16,
     7 8.17900786477396d-16/
      data db(1),  db(2),  db(3),  db(4),  db(5),  db(6),  db(7),
     1     db(8),  db(9),  db(10), db(11), db(12), db(13), db(14),
     2     db(15)          /-2.77571356944231d-01, 4.44212833419920d-03,
     3-8.42328522190089d-05,-2.58040318418710d-06, 3.42389720217621d-07,
     4-6.24286894709776d-09,-2.36377836844577d-09, 3.16991042656673d-10,
     5-4.40995691658191d-12,-5.18674221093575d-12, 9.64874015137022d-13,
     6-4.90190576608710d-14,-1.77253430678112d-14, 5.55950610442662d-15,
     7-7.11793337579530d-16/
c***first executable statement  djairy
      if (x.lt.0.0d0) go to 90
      if (c.gt.5.0d0) go to 60
      if (x.gt.1.20d0) go to 30
      t = (x+x-1.2d0)*con4
      tt = t + t
      j = n1
      f1 = ak1(j)
      f2 = 0.0d0
      do 10 i=1,m1
        j = j - 1
        temp1 = f1
        f1 = tt*f1 - f2 + ak1(j)
        f2 = temp1
   10 continue
      ai = t*f1 - f2 + ak1(1)
c
      j = n1d
      f1 = dak1(j)
      f2 = 0.0d0
      do 20 i=1,m1d
        j = j - 1
        temp1 = f1
        f1 = tt*f1 - f2 + dak1(j)
        f2 = temp1
   20 continue
      dai = -(t*f1-f2+dak1(1))
      return
c
   30 continue
      t = (x+x-con2)*con3
      tt = t + t
      j = n2
      f1 = ak2(j)
      f2 = 0.0d0
      do 40 i=1,m2
        j = j - 1
        temp1 = f1
        f1 = tt*f1 - f2 + ak2(j)
        f2 = temp1
   40 continue
      rtrx = sqrt(rx)
      ec = exp(-c)
      ai = ec*(t*f1-f2+ak2(1))/rtrx
      j = n2d
      f1 = dak2(j)
      f2 = 0.0d0
      do 50 i=1,m2d
        j = j - 1
        temp1 = f1
        f1 = tt*f1 - f2 + dak2(j)
        f2 = temp1
   50 continue
      dai = -ec*(t*f1-f2+dak2(1))*rtrx
      return
c
   60 continue
      t = 10.0d0/c - 1.0d0
      tt = t + t
      j = n1
      f1 = ak3(j)
      f2 = 0.0d0
      do 70 i=1,m1
        j = j - 1
        temp1 = f1
        f1 = tt*f1 - f2 + ak3(j)
        f2 = temp1
   70 continue
      rtrx = sqrt(rx)
      ec = exp(-c)
      ai = ec*(t*f1-f2+ak3(1))/rtrx
      j = n1d
      f1 = dak3(j)
      f2 = 0.0d0
      do 80 i=1,m1d
        j = j - 1
        temp1 = f1
        f1 = tt*f1 - f2 + dak3(j)
        f2 = temp1
   80 continue
      dai = -rtrx*ec*(t*f1-f2+dak3(1))
      return
c
   90 continue
      if (c.gt.5.0d0) go to 120
      t = 0.4d0*c - 1.0d0
      tt = t + t
      j = n3
      f1 = ajp(j)
      e1 = ajn(j)
      f2 = 0.0d0
      e2 = 0.0d0
      do 100 i=1,m3
        j = j - 1
        temp1 = f1
        temp2 = e1
        f1 = tt*f1 - f2 + ajp(j)
        e1 = tt*e1 - e2 + ajn(j)
        f2 = temp1
        e2 = temp2
  100 continue
      ai = (t*e1-e2+ajn(1)) - x*(t*f1-f2+ajp(1))
      j = n3d
      f1 = dajp(j)
      e1 = dajn(j)
      f2 = 0.0d0
      e2 = 0.0d0
      do 110 i=1,m3d
        j = j - 1
        temp1 = f1
        temp2 = e1
        f1 = tt*f1 - f2 + dajp(j)
        e1 = tt*e1 - e2 + dajn(j)
        f2 = temp1
        e2 = temp2
  110 continue
      dai = x*x*(t*f1-f2+dajp(1)) + (t*e1-e2+dajn(1))
      return
c
  120 continue
      t = 10.0d0/c - 1.0d0
      tt = t + t
      j = n4
      f1 = a(j)
      e1 = b(j)
      f2 = 0.0d0
      e2 = 0.0d0
      do 130 i=1,m4
        j = j - 1
        temp1 = f1
        temp2 = e1
        f1 = tt*f1 - f2 + a(j)
        e1 = tt*e1 - e2 + b(j)
        f2 = temp1
        e2 = temp2
  130 continue
      temp1 = t*f1 - f2 + a(1)
      temp2 = t*e1 - e2 + b(1)
      rtrx = sqrt(rx)
      cv = c - fpi12
      ccv = cos(cv)
      scv = sin(cv)
      ai = (temp1*ccv-temp2*scv)/rtrx
      j = n4d
      f1 = da(j)
      e1 = db(j)
      f2 = 0.0d0
      e2 = 0.0d0
      do 140 i=1,m4d
        j = j - 1
        temp1 = f1
        temp2 = e1
        f1 = tt*f1 - f2 + da(j)
        e1 = tt*e1 - e2 + db(j)
        f2 = temp1
        e2 = temp2
  140 continue
      temp1 = t*f1 - f2 + da(1)
      temp2 = t*e1 - e2 + db(1)
      e1 = ccv*con5 + 0.5d0*scv
      e2 = scv*con5 - 0.5d0*ccv
      dai = (temp1*e1-temp2*e2)*rtrx
      return
      end