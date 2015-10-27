*deck @(#)iounq.f	5.1  11/6/94
      character*(*) function iounq()
c
c***begin prologue     iounq
c***date written       860112   (yymmdd)
c***revision date      860112   (yymmdd)
c***keywords           iosys dependent routines
c
c***author             saxe, paul,    (lanl)
c***source             @(#)iounq.f	5.1   11/6/94
c***purpose            to find a physical file name for a scratch file.
c                         the file must not already exist.
c***description        #
c
c
c***references
c
c***routines called    
c
c   common blocks:     (none)
c
c***end prologue       iounq
c
      implicit integer (a-z)
c
      character char*5
      character*8 name
      logical ioinq
c
      data char /'vwxyz'/
c
c
      name='t/'
      do 3 i=1,5
         name(5:5)=char(i:i)
         do 2 j=1,5
            name(6:6)=char(j:j)
            do 1 k=1,5
               name(7:7)=char(k:k)
               if (.not.ioinq(name,'exist')) then
                  iounq=name
                  return
               end if
    1       continue
    2    continue
    3 continue
c
      call lnkerr('cannot create a unique filename!!!')
c
      return
      end
