*deck @(#)sorter.f	5.1  11/6/94
      subroutine sorter(entry,icore,core,maxt,n,labels,bins,values,
     #                  npasst,filet,unitt,prnt)
c
c***begin prologue     sorter
c***date written       850601   (yymmdd)
c***revision date      870216   (yymmdd)
c
c   16 february 1987   pws at lanl
c       for the one-pass sort, the iosys call to create 'bins' was missing
c       the last argument.
c
c***keywords           sorting, yoshimine sort
c
c***author             saxe, paul,    (lanl)
c***source             @(#)sorter.f	5.1   11/6/94
c***purpose            to sort, on disc, or in core, a real numbers
c                        with an associted desired position.
c***description        sorter is intended to receive as input real
c numbers and an associated position from 1 to n; however, this input
c is in random order. finally, sorter will place a file containing
c the sorted real numbers (according to the associated tag) on a
c specified iosys unit. the call has the following form:
c
c    call sorter(directive,icore,core,i1,i2,i3,i4,r1,i5,c1,c2)
c
c the arguments i1-i5 are integer arguments, r1 is real, and c1 and c2
c are charcter arguments. the actual use of these arguments depends on
c the directive, and will be described in detail later. 'directive' is
c a character string which directs tha action sorter will take. 'icore'
c and 'core' are integer and real scratch areas, respectively, which
c are assumed to occupy the same physical memory. i.e. 'icore' and
c 'core' are equivalenced to each other either explicitly or
c implicitly.
c
c      the possible directives are:
c              'start'
c              'with bin'
c              'end'
c
c each of these is now described in detail:
c
c call sorter('start',icore,core,maxcor,n,minbin,accume,r1,npass,
c              file,unit)
c
c where:
c       maxcor   integer
c                the length of icore in integer words.
c                the length of core is iadtwp(maxcor) real words
c
c       n        integer
c                the maximum number of data to sort, i.e. the maximum
c                value of the tag, or final location, of any datum. the
c                data is assumed to occupy real*8 words.
c
c       minbin   integer
c                the minimum bin size allowable in the sort. defaults
c                to 512 elements.
c
c       accume   integer
c                .ge.0 to zero the final array and put the elements
c                      in their places, so that elements not present
c                      in the sort input will be zero.
c                .lt.0 to read in the final array from the final file,
c                      and add the sorted values in to their locations
c                      current value.
c
c       r1       real array
c                unreferenced.
c
c       npass    integer
c                if non-zero, a value of the number of passes to force
c                in the sort. currently possible values are
c                      -1  in-core sort
c                       1  one pass sort
c                       2  two pass sort
c
c       file     character
c                the name of the iosys file where the final array will
c                reside, and also, if accume=-1, where the initial
c                array is found.
c
c       unit     character
c                the name of the iosys unit for the final file.
c
c        this call is used to initialize the various paramters for the
c   sort. at this point, sorter decides how to do the sort in the space
c   and restrictions given. if sorter decides that is unable to do the
c   sort in the space given, it ungraciously pulls the plug by calling
c   lnkerr. if, however, the sort is feasible, sorter initializes its
c   internal arrays, etc.
c
c
c call sorter('with bin',icore,core,i1,nvals,labels,bins,values,i2,c1,
c              c2)
c
c where:
c       nvals    integer
c                the number of data being passed in in this call.
c
c       labels   integer (nvals)
c                the tags for the data, i.e. the final location for
c                each datum.
c
c       bins     integer (nvals)
c                a scratch array, unreferenced in an in-core sort.
c
c       values   real (nvals)
c                the data to be sorted.
c
c       i1,i2    integer
c                unreferenced.
c
c       c1,c2    character
c                unreferenced.
c
c        this call passes in data and the associated tags to sorter. on
c   a vector machine it helps to keep 'nvals' reasonably long.
c
c
c call sorter('end',icore,core,i1,i2,i2,i4,i5,r1,i6,c1,c2)
c
c         this call instructs sorter to complete the sort of the data
c     submitted 'with bin'. the sorted data is written to an iosys
c     file described by 'file' and 'unit' in the initialization call.
c     the arguments i1-i6,r1,c1 and c2 are unreferenced integer, real
c     and character arguments. i4, i5 and r1 are arrays.
c
c***references
c
c***routines called    lnkerr (mdutil)
c                      rzero  (math)
c                      iosys  (io)
c                      unqfil (math)
c                      sorter1 (math)
c                      sorter2 (math)
c                      sorter3 (math)
c                      sorter4 (math)
c                      scatter (calmath)
c
c   common blocks:     io
c
c   iosys units:       sort   (unqfil('sort'))  1 and 2 pass sorts only
c                      sort1  (unqfil('sorta')) 2 pass sorts only
c
c***end prologue       sorter
c
      implicit integer (a-z)
c
      character*(*) entry,filet,unitt
      character*32 file
      character*16 unit
      character*8 cjunk
      character*1 itoc
      real*8 core(*),values(*)
      integer icore(*),labels(*),bins(*)
      logical prnt,debug
c      logical writit
c
      parameter (debug=.true.)
c
      common /io/     inp,iout
c      common/test/writit
c
      save ndata,minbin,accume,file,unit,maxcor,npass,nbins,binsiz
      save chain,lbin,vbin,perbin,offset,bintot,nsortd
      save perbn1,perbn2,binsz1,binsz2
c
c
c
      if (entry.eq.'start') then
         ndata=n
         minbin=(max(512,labels(1))+511)/512*512
         accume=bins(1)
         file=filet
         unit=unitt
         maxcor=maxt
         offset=0
         bintot=0
         nsortd=0
c
c        ----- determine how many passes will be needed -----
c
         if (maxcor.ge.wptoin(ndata)) then
            npass=0
         else
            perbin=(maxcor-(1+wptoin(1))*minbin)/512*512+512
            do 1 junk=1,100
               nbins=(ndata-1)/perbin+1
               perbin=(maxcor-nbins-(1+wptoin(1))*minbin)/512*512+512
               if ((ndata-1)/perbin+1.eq.nbins) go to 2
    1       continue
c
            call lnkerr('internal problem in sorter with nbins')
c
    2       continue
            npass=(nbins+(1+wptoin(1))*minbin*nbins)/maxcor+1
         end if
         if (debug) then
            write(iout,*) 'Initial sort information'
            write(iout,*)
            write(iout,*) 'ndata   = ',ndata,'maxcor  = ',maxcor,
     $                    'nbin    = ',nbins
            write(iout,*)
            write(iout,*) 'perbin  = ',perbin,'npass    = ',npass
         endif
c
c        ----- check feasibility of sort -----
c
         if (npass.gt.2) call lnkerr('sorter cannot do sort in less'//
     #                               ' than three passes')
c
c        ----- check for overide of number of passes -----
c
         if (npasst.ne.0) then
            if (npasst.eq.-1) then
               if (npass.ne.0) call lnkerr('cannot force in core sort')
               npass=0
            else
               if (npasst.lt.npass) call lnkerr('cannot force '//
     #                 itoc(npasst)//' pass sort')
               npass=npasst
            end if
         end if
c
c        ----- allocate core for bins, etc. -----
c
         chain=1
c
         if (npass.eq.0) then
c
c           ----- in-core sort is a special case, and will be handled
c                 separately
c
c           ----- rewind or create output file -----
c
            if (accume.ge.0) then
                call iosys('get maximum length of "'//file//'" on '//
     #                      unit,len,0,0,' ')
               if (len.le.0) then
                  call iosys('create real "'//file//'" on '//unit,
     #                       ndata,0,0,' ')
               else
                  if (len.lt.ndata) then
                     call lnkerr(file//' on '//unit//' already '//
     #                           'exists but is not long enough')
                  end if
               end if
               call rzero(core,ndata)
            else
               call iosys('read real "'//file//'" from '//unit,ndata,
     #                     core,0,' ')
            end if
c
            if (prnt) then
               write (iout,3)
    3          format (5x,'sorter will use in-core algorithm')
            end if
c
         else if (npass.eq.1) then
            if(debug) then
               write(iout,*) 'One pass sort: checking max. bin size'
            endif
c
c            ----- find the largest possible bin size -----
c
            do 4 i=2,max(2,nbins+24)
               binsiz=((maxcor-i)/(i+wptoin(i)))/512*512
               perbin=(ndata/i+1)/512*512+512
               size1=binsiz+wptoin(binsiz)+wptoin(perbin)+i
               size2=i*(binsiz+wptoin(binsiz))+i
               if (size1.le.maxcor.and.size2.le.maxcor) go to 5
    4       continue
c
            call lnkerr('internal failure in sorter for 1 pass')
c
    5       if(debug) then
               write(iout,*) 'One pass sort information'
               write(iout,*)
               write(iout,*) 'binsiz  = ',binsiz,'perbin  = ',perbin  
               write(iout,*)
               write(iout,*) 'size1   = ',size1,'size2   = ',size2  
            endif
            nbins=i
            lbin=chain+nbins
            vbin=iadtwp(lbin+binsiz*nbins)
c
c           ----- initialize the bins -----
c
            do 6 i=1,nbins
               icore(lbin+(i-1)*binsiz  )=2
               icore(lbin+(i-1)*binsiz+1)=-1
    6       continue
c
c           ----- and open a temporary sort file -----
c
c.ltss            len=max(32768,ndata/20)
c.ltss            len=min(len,16777216)
c
            len=2*ndata+200000
            call iosys('open sort as scratch on ssd',len,0,0,cjunk)
            call iosys('create integer bins on sort',-1,0,0,' ')
c
            if (prnt) then
               write (iout,7) nbins,binsiz
    7          format (5x,'sorter will use one pass of ',i4,' bins of ',
     #                   i8,' words each')
            end if
c
c
         else if (npass.eq.2) then
c
c           ----- find maximum bin size -----
c
            do 8 i=2,max(2,nbins+10)
               binsz1=((maxcor-i)/(i+wptoin(i)))/512*512
               perbn1=ndata/i+1
               binsz2=((maxcor-i-binsz1-wptoin(binsz1))/(i+wptoin(i)))/
     #                 512*512
               perbn2=(perbn1/i+1)/512*512+512
               size1=binsz2+wptoin(binsz2)+wptoin(perbn2)+2*i
               size2=i*(binsz1+wptoin(binsz1))+2*i
               size3=i*(binsz2+wptoin(binsz2))+2*i+2*binsz1+
     #               wptoin(binsz1)
               if (size1.le.maxcor.and.size2.le.maxcor.and.
     #             size3.le.maxcor) go to 9
    8       continue
c
            call lnkerr('internal failure in sorter for 2 pass sort')
c
    9       continue
c
c           ----- allocate core -----
c
            nbins=i
            lbin=chain+nbins
            vbin=iadtwp(lbin+binsz1*nbins)
            top=wpadti(vbin+binsz1*nbins)
c
            if (top.gt.maxcor) call lnkerr('internal error in sorter'//
     #          ' allocating core for first phase of two pass sort')
c
c           ----- initialize the bins -----
c
            do 10 i=1,nbins
               icore(lbin+(i-1)*binsz1  )=2
               icore(lbin+(i-1)*binsz1+1)=-1
   10       continue
c
c           ----- and open temporary sort files -----
c
c.ltss            len=max(32768,perbn1/20)
c.ltss            len=min(len,16777216)
c
            len=2*perbn1+200000
            call iosys('open sort1 as scratch on ssd',len,0,0,cjunk)
            call iosys('create integer bins on sort1',-1,0,0,' ')
c
c           ----- try to guarantee enough space on the ssd -----
c
            do 11 i=1,(perbn1+wptoin(perbn1))/maxcor+1
               len=min(perbn1+wptoin(perbn1),maxcor)
               call iosys('write integer bins on sort1 '//
     #                    'without rewinding',len, core,0,' ')
   11       continue
            call iosys('rewind bins on sort1',0,0,0,' ')
c
c           ----- and create first phase sort tape -----
c
c.ltss            len=max(32768,ndata/20)
c.ltss            len=min(len,16777216)
c
            len=2*ndata+200000
            call iosys('open sort as scratch on ssd',len,0,0,cjunk)
            call iosys('create integer bins on sort',-1,0,0,' ')
c
            if (prnt) then
               write (iout,12) nbins,binsz1,perbn1,binsz2,perbn2
   12          format (5x,'sorter will use three passes of ',i4,' bins'
     #                /,7x,i8,' words each and ',i8,' elements ',
     #                'for phase 1',/,7x,i8,' words each and ',i8,
     #                ' elements for phase 2')
            end if
c
c           ----- transfer sizes to general variables -----
c
            binsiz=binsz1
            perbin=perbn1
c
         else
            call lnkerr('internal error in sorter with core')
         end if
         return
c
c     ----- entry with values to sort -----
c
      else if (entry.eq.'with bin') then
c
         if (npass.eq.0) then
c
c           ----- incore sort -----
c
            if (accume.ge.0) then
               call scatter(n,core,labels,values)
            else
               do 13 i=1,n
                  core(labels(i))=core(labels(i))+values(i)
   13          continue
            end if
         else
c
c           ----- out-of-core sort -----
c
            call sorter1(n,labels,bins,values,icore(lbin),core(vbin),
     #                   binsiz,perbin,offset,bintot,'sort')
         end if
         nsortd=nsortd+n
c
         return
c
c     ----- finish the sorting -----
c
      else if (entry.eq.'end') then
c         write(iout,*) 'The entry variable is', entry
c         write(iout,*) 'npass ndata', npass, ndata
c
         if (npass.eq.0) then
c            write(iout,*) 'Doing an incore'
c
c           ----- in-core sort -----
c
            call iosys('write real "'//file//'" on '//unit,ndata,
     #                  core,0,' ')
c
         else if (npass.eq.1) then
c
c           ----- flush the final partial bins -----
c
c            write(iout,*) 'using sorter2 on a 1 pass'
c            writit=.true.
            call sorter2(icore(lbin),core(vbin),binsiz,nbins,
     #                   core(chain),bintot,'sort')
c            writit=.false.
c
c           ----- write eof on bins -----
c
            call iosys('endfile bins on sort',0,0,0,' ')
c
            if (prnt) then
               write (iout,14) nsortd,bintot
   14          format (5x,'one pass sort wrote ',i9,' elements in ',i6,
     #                 ' bins')
            end if
c
c           ----- reallocate core and loop back through bins -----
c
            lbin=chain+nbins
            vbin=iadtwp(lbin+binsiz)
            sorted=vbin+binsiz
            top=wpadti(sorted+perbin)
c
            if (top.gt.maxcor) then
               call lnkerr('internal error on second phase of one'//
     #                     ' pass sort')
            end if
c
            call sorter3(nbins,binsiz,icore(chain),icore(lbin),
     #                   core(vbin),core(sorted),perbin,accume,
     #                   file,unit,ndata,'sort',offset)
c
            call iosys('destroy sort',0,0,0,' ')
c
c
            if (prnt) then
               write (iout,15) ndata,file,unit
   15          format (5x,'one pass sort wrote ',i9,' elements to ',a16,
     #                   ' on ',a16)
            end if
c
         else if (npass.eq.2) then
c            write(iout,*) 'Its 2 passes'
c
c           ----- flush the final partial bins -----
c
            call sorter2(icore(lbin),core(vbin),binsiz,nbins,
     #                   core(chain),bintot,'sort')
c
c           ----- write eof on bins -----
c
            call iosys('endfile bins on sort',0,0,0,' ')
c
            write (iout,16) nsortd,bintot
   16       format (5x,'first pass wrote ',i9,' elements in ',i6,
     $           ' bins')
c
c           ---- reallocate core -----
c
            chain2=chain+nbins
            lbin=chain2+nbins
            labs=lbin+binsz2*nbins
            bin=labs+binsz1
            vbin=iadtwp(bin+binsz1)
            vals=vbin+binsz2*nbins
            top=wpadti(vals+binsz1)
c
            lab2=chain2+nbins
            val2=iadtwp(lab2+binsz2)
            sorted=val2+binsz2
            top2=wpadti(sorted+perbn2)
c
            if (top.gt.maxcor.or.top2.gt.maxcor) then
               call lnkerr('internal core allocation error for '//
     #                     'second pass of two-pass sort')
            end if
c
c
            call sorter4(core(vals),core(vbin),core(val2),core(sorted),
     #                   icore(labs),icore(lbin),icore(bin),
     #                   icore(lab2),icore(chain),icore(chain2),
     #                   binsz1,binsz2,nbins,perbn1,perbn2,accume,
     #                   file,unit,ndata)
c
            call iosys('destroy sort',0,0,0,' ')
            call iosys('destroy sort1',0,0,0,' ')
c
            if (prnt) then
               write (iout,17) ndata,file,unit
   17          format (5x,'two pass sort wrote ',i9,' elements to ',a16,
     #                   ' on ',a16)
            end if
c
         else
            call lnkerr('bad number of passes as finish sort')
         end if
c
c
         return
c
c     ----- handle strange entries -----
c
      else
         unit=entry
         call lnkerr('error with entry to sorter:'//unit)
      end if
c
      end
