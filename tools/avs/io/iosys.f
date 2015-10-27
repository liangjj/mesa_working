*deck @(#)iosys.f	4.1  7/7/93
      subroutine iosys(str,i1,array,i2,c1)
c
c***begin prologue     iosys
c***date written       850125   (yymmdd)
c***revision date      870201   (yymmdd)
c
c   1 february 1987  pws at lanl
c      ioopen: changing 'unitnm' and 'number'to character*256 variables to
c      allow for long file names to accomodate unix directory system.
c
c   31 january 1987  pws at lanl
c
c      bsd4.2 unix version on sun microsystems 3/50 and 3/160.
c      had to modify the fortran open statements in ioopen to 
c      remove some vaxism's such as 'initialsize'. also changed
c      'idum' from integer to character*8 and made 'pt' and 'nfiles'
c      character when putting them into or getting them from 'idum'.
c      finally, modified ioexist ('does ... exist on ...') to return
c      a character answer, rather than hollerith, since f77 complains
c      about an integer being set equal to a character string.
c
c***keywords           io, i/o, files, input, output
c
c***author             saxe, paul,    (lanl)
c***source             @(#)iosys.f	4.1   7/7/93
c***purpose            flexible random-access unformatted i/o.
c
c***description        iosys is an attempt at a flexible, user-friendly
c i/o package which hides the machine dependent aspects of i/o. all i/o
c in the mesa system other than fortran formatted input and output
c should be performed using iosys. the form of a call to iosys is:
c
c          call iosys(directive,arg1,arg2,arg3,arg4)
c
c where 'directive' is a character string giving (in english?) the
c action to be taken, and the remaining three arguments depend on what
c is being done.
c      in the following description of the directives <file> is an
c iosys internal file name; <unit> is an iosys internal unit name;
c parenthesis () indicate a choice of words, one of which must be used;
c and brackets [] indicate an optional keyword or phrase. also, all
c (hopefully) non-essential words such as 'to', 'on' and 'of' may be
c omitted wherever they appear; however, their inclusion is suggested
c as an aid to readability.
c      the possible directives are:
c        'length of <file> on <unit>'
c        'get length of <file> on <unit>'
c        'get maximum length of <file> on <unit>'
c        'get read pointer of <file> on <unit>'
c        'get write pointer of <file> on <unit>'
c        'does <file> exist on <unit>'
c        'write (character integer real) <file> (to,on) <unit>
c                    [without rewinding] [asynchronously]'
c        'read (character integer real) <file> from <unit>
c              [without rewinding] [asynchronously]'
c        'end <file> on <unit>'
c        'endfile <file> on <unit>'
c        'create (character integer real) <file> on <unit>'
c        'open <unit> [as (new,old,unknown,scratch)] [on ssd]'
c        'dump'
c        'destroy [<unit>]'
c        'close (all,<unit>)'
c        'wait for (all,<unit>)'
c        'rewind (all,<file>) on (all,<unit>)
c                 [(read,write,read-and-write)]'
c        'set alignment (for,of) <unit>'
c        'reset alignment (for,of) <unit>'
c
c      iosys manages collections of 'files' on 'units', which are
c actually physical disc files. iosys uses a directory to keep track
c of its logical 'files' on each of its logical 'units'. the first step
c in using iosys is to open one or more logical units using the 'open'
c directive. this is the only stage where the user must be cognizant of
c the physical file name associated with the logical units. after a
c unit is opened, it is always referred to by its logical name
c embedded in a directive string.
c
c     there are two ways to create a logical file on a unit: the first
c is simply by writing data to a new file name using the 'write'
c directive, the second is with an explicit 'create' directive. the
c first method will create a file exactly the length of the data being
c written; whereas, a 'create' directive allows the user to specify
c the length of the file. as a user, you will only need to explicitly
c create files only if all the data cannot be written to the file in
c one write. there is also a variant of the explicit create invoked
c by giving a length of -1. in this case, an 'endless' file is created,
c on which an unspecified amount of data may be written. an 'endless'
c file is active until an 'end' or 'endfile' directive is issued, at
c which point the 'endless' file is fixed to the current length and may
c not be further extended. there are restrictions on the use of
c 'endless' files: first, only one endless file can be active on a
c given unit at one time; and second, no new files may be created
c implicitly or explicitly on a unit with an active endless file.
c
c      logical file and unit names may have up to 16 characters; any
c beyond 16 are truncated. if the name contains embedded blanks, it
c must be enclosed in double quotes, e.g. "my file". the quotes are
c treated as part of the name, so a name surrounded with quotes is
c not the same as the name without quotes. iosys can currently handle
c 20 units open simultaneously with a total of 1000 files. these limits
c are set in the parameter statement and are readily modified.
c
c      the length of a file or of a transfer (read or write) are in
c terms of the type of data on the file. integer and real files are
c in terms of integer or real words respectively, while character
c data is by character. also note that for character data, the length
c need not be given as it can be assumed to be the length of the
c character string passed in.
c
c      'alignment' has to do with the physical storage of data. if
c alignment is set, new files created will be positioned starting at
c physical disc sectors. this will speed up i/o operations, but will
c also make the physical file larger since gaps will be left between
c logical files. another ramification of alignment is that the
c 'maximum length' of files can be larger than their 'length' since
c iosys will allow a file to expand in to the gap left by alignment.
c alignment is turned off with a 'reset alignment' directive. by
c default, units are not aligned. finally, the alignment of a unit can
c be changed as often as desired since it only affects new files being
c created.
c
c      this section will describe briefly each directive and the other
c arguments associated with it. a 'u' will denote that the argument is
c not referenced when using a directive, and the user should use a '0'
c or any variable in the calling sequence.
c
c manipulating units:
c
c call iosys('open <unit> [as (new,old,unknown)] [on ssd]',size,u,
c            u,name)
c
c           'name' is the hollerith (or character) name or a unit
c           number of the physical disc file to use. it must have 7
c           or less characters.
c
c           'size' is the size of family members on ctss systems.
c
c            the default for an open is 'unknown', which means that if
c            the file exists, it will be opened, if not, it will be
c            created. 'old' specifies that the file must exist;
c            otherwise, there is an error. 'new' guarantees a new file;
c            if it already exists, it will be deleted and then
c            recreated. 'scratch' can be used to create a unique
c            temporary file. in this case, no physical file name is
c            needed as iosys will create one, and when the unit is
c            closed, the physical file will be deleted.
c            'on ssd' specifies that the file should be put
c            on a solid-state disc if possible. if there is not room,
c            or no ssd, the file is created on a normal disc. the ssd
c            should only be used for heavily used, temporary files
c            that only exist as scratch in one program link.
c
c
c call iosys('close (all,<unit>)',u,u,u,u)
c
c           this directive closes either a specific unit or all open
c           units. when a unit is closed, the directory information
c           is saved on the unit so that a subsequent open will restore
c           the unit to the identical state as when it was closed--
c           including the read and write pointers for files.
c
c
c call iosys('dump',u,u,u,u)
c
c           this will dump a listing of the file names, pointers, etc
c           for all units to the output file.
c
c call iosys('destroy [<unit>]',u,u,u,name)
c
c           this directive works in two fashions: if an iosys unit name
c           is given in the directive, that iosys unit is closed and
c           the physical file is destroyed. in this case, the argument
c           'name' is not referenced. otherwise, the physical file
c           'name' is destroyed, regardless of whether it is an iosys
c           file or not.
c
c call iosys('set alignment on <unit>',u,u,u,u)
c call iosys('reset alignment on <unit>',u,u,u,u)
c
c           these calls turn on and off the alignment of new files to
c           physical sector boundaries.
c
c
c manipulating files:
c
c call iosys('create (character integer real ) <file> on <unit>',
c            length,u,u,u)
c
c           this call creates a file of length 'length' data elements
c           on an iosys unit. if 'length' is -1, an 'endless' file
c           is created (see above for a description of endless files)
c
c call iosys('write (character integer real) <file> (to,on) <unit>
c            [without rewinding] [asynchronously]',nwords,array,offset,
c            c)
c
c           'nwords' is the number of data elements to transfer.
c           'array' is the array to transfer from.
c           'offset' is added to the current write pointer before
c               writing.
c           'c' is the character string if writing characters.
c
c           if 'nwords' is -1, iosys will use the current size of the
c           file and transfer that number of data elements, which
c           is a nice shortcut for rewriting a file. the default is
c           to rewind the file before writing. if this is not overidden
c           with 'without rewinding', then 'offset' specifies a
c           location in the file to begin writing (starting at 0).
c           otherwise, 'offset' is added to the current write location
c           before writing. thus, repeated calls 'without rewinding'
c           and an 'offset' of 0 will sequentially write a file. note
c           that the read and write pointers are separate and do not
c           affect each other. for character data, the length need
c           not be given, in which case the length of 'c' will be used.
c
c call iosys('read (character integer real) <file> from <unit>
c            [without rewinding] [asynchronously]',nwords,array,offset,
c            c)
c
c           reading is exactly the inverse of writing, and all options
c           and arguments are identical to the above.
c
c call iosys('end <file> on <unit>',u,u,u,u)
c call iosys('endfile <file> on <unit>,u,u,u,u)
c
c           use this to terminate an 'endless' file. the length of the
c           file is fixed at the length at the time of this call.
c
c call iosys('wait for <unit>',u,u,u,u)
c
c           used to wait for the completion of an asynchronous i/o
c           request. an array being written asynchronously cannot
c           safely be modified until a 'wait' call has been made.
c           converseley, data being read asynchronously cannot be used
c           before a call to 'wait'.
c
c call iosys('rewind (all,<file>) on (all,<unit>) [(read,write,
c            read-and-write)]',u,u,u,u)
c
c           rewinds the read and/or write pointers of one or all of
c           the files on one or all units.
c
c
c inquiries about files:
c
c call iosys('length of <file> on <unit>',length,u,u,u)
c call iosys('get length of <file> on <unit>',length,u,u,u)
c
c           returns the length of the data written
c           to a file. if the file does not exist, -1 is returned.
c
c call iosys('get maximum length of <file> on <unit>',length,u,u,u)
c
c           the maximum amount of data which can be
c           stored on a file. if the file does not exist, a -1 is
c           returned.
c
c call iosys('get read pointer of <file> on <unit>',pointer,u,u,u)
c call iosys('get write pointer of <file> on <unit>',pointer,u,u,u)
c
c           these are the addresses in a file (counting
c           from 0) where the next read or write would occur if the
c           file is not rewound. if the file does not exist, a -1 is
c           returned.
c
c call iosys('does <file> exist on <unit>',u,u,u,answer)
c
c           'answer' returns as 'no' or 'yes'.
c
c
c
c
c
c***references
c
c***routines called    gettok (io)
c                      ioopen (io)
c                      iodest (io)
c                      iodump (io)
c                      iofile (io)
c                      ioeof  (io)
c                      iowrit (io)
c                      ioread (io)
c                      ioclos (io)
c                      iowait (io)
c                      iorew  (io)
c                      iounit (io)
c                      lnkerr (mdutil)
c                      ioexst (io)
c                      iolen  (io)
c
c   common blocks:     ioqqq1, ioqqq2, and io in iodump.
c
c***end prologue       iosys
c
      implicit integer (a-z)
c
      character*(*) str,c1
c
      character*16 token
      integer array(1)
c
      character*240 string
      character*80 tmplin
      character*32 keylis,key
      character*16 cconst,unlist,unityp,type,unit
      character*8  filtyp
      integer base,readpt,writpt,eof,end,iconst,unitpt,nfile
      integer un,file,nunits,unitno
      real*8 rconst
      logical align,locked,called
c
c     ----- functions -----
c
      character*240 dcaptl
c
      parameter (maxfil=2000,maxun=20)
c
      common /ioqqq1/ keylis(maxfil),filtyp(maxfil),cconst(maxfil),
     #                unlist(maxun),type,key,unit,tmplin,unityp(maxun)
      common /ioqqq2/ base(maxfil),readpt(maxfil),writpt(maxfil),
     #                eof(maxfil),end(maxfil),iconst(maxfil),
     #                rconst(maxfil),unitpt(maxun),nfile(maxun),
     #                align(maxun),locked(maxun),unitno(maxun),
     #                un,file,nunits
c
c
      data called/.false./
      save called
c
c     set up some common block information.
      if(.not.called) then
         nunits=0
         do 10 i=1,maxun
            align(i)=.true.
            locked(i)=.false.
   10    continue
         called=.true.
      end if
c
c     ----- transfer the command string to tmplin in case we need
c           to print the string for an error
c
      tmplin=str
c
c     ----- convert the command string to lower case -----
c
      string=dcaptl(str)
c
c     ----- get the first token from 'string' -- the operation -----
c
      pos=0
      call gettok(token,string,pos)
c
c     ----- decide what to do based on 'oper' -----
c
      if (token.eq.'read') then
         call ioread(string,pos,i1,array,array,i2,c1)
      else if (token.eq.'write') then
         call iowrit(string,pos,i1,array,array,i2,c1)
      else if (token.eq.'copy') then
         call iocopy(string,pos,array,i1,c1)
      else if (token.eq.'open') then
         call ioopen(string,pos,c1,i1)
      else if (token.eq.'destroy') then
         call iodest(string,pos,i1)
      else if (token.eq.'dump') then
         call iodump
      else if (token.eq.'create') then
         call iofile(string,pos,i1)
      else if (token.eq.'end'.or.token.eq.'endfile') then
         call ioeof(string,pos)
      else if (token.eq.'close') then
         call ioclos(string,pos)
      else if (token.eq.'wait') then
         call iowait(string,pos)
      else if (token.eq.'rewind') then
         call iorew(string,pos)
      else if (token.eq.'set') then
         call gettok(token,string,pos)
         if (token(1:5).eq.'align') then
            call gettok(token,string,pos)
            if (token.eq.'for'.or.token.eq.'of') call gettok(token,
     #                                           string,pos)
            un=iounit(token,unlist,nunits,error)
            if (error.ne.0) then
               call lnkerr('io system: trying to set alignment on'
     #                     //' non-existant unit '//token)
         end if
            align(un)=.true.
         else
            call lnkerr('io system: do not understand what to set: '
     #                  //tmplin)
         end if
      else if (token.eq.'reset') then
         call gettok(token,string,pos)
         if (token(1:5).eq.'align') then
            call gettok(token,string,pos)
            if (token.eq.'for'.or.token.eq.'of') call gettok(token,
     #                                           string,pos)
            un=iounit(token,unlist,nunits,error)
            if (error.ne.0) then
               call lnkerr('io system: trying to reset alignment on'
     #                     //' non-existant unit '//token)
            end if
         else
            call lnkerr('io system: do not understand what to reset:'
     #                  //tmplin)
         end if
      else if (token.eq.'does') then
         call ioexst(string,pos,c1)
      else if (token.eq.'length'.or.token.eq.'get') then
         call iolen(string,pos,i1)
      else if (token.eq.'delete') then
      else
         call lnkerr('io system: unrecognised operation: '//
     #                tmplin)
      end if
c
      return
      end
