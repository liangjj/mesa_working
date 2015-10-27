2# Repository for Mesa sources, libraries and scripts

     * All of the mesa links are in this repository.
     * All of the library files needed by mesa as well as many subroutines which have been written over 30 years reside here.
     * Scripts needed to set up the environment and define machine dependent quantities are also here.
     * There are a number of versions of mesa that exist.  The version in this repository has been continuously updated
       for many years.  The library is far more extensive than the original library as far as the number of routines
       that are available to the user. In the original source code there was a provision for dynamic allocation of
       arrays but that was never really implemented.  A large blank common was used.  This has been replaced in a
       number of the links by true dynamic allocation.  In addition, the memory requested in those links was carefully
       looked at.  In some of the links, the original code was quite sloppy and that needs more work.  Another major
       change was the complete re-write of the Davidson algorithm to avoid the disk thrashing that was going on in the
       original code.  Now memory is allocated for all vectors and iterates.  The routine is far more efficient than
       before. Since the code has not been used in a long time it does need to be tested.

## Installing MESA

    1.  Clone the scripts_linux and/or the scripts_mac folders from Bitbucket.  Put them in your home directory.
    2.  In the scripts folder there is a copy of a .bash_profile.  Open that file with your favorite editor.  You will see that various
        variables are defined and exported.  Two are important; the location of the scripts folder and the location of the
        working_mesa_directory.  I have already addressed the former.  My strong suggestion is to create a folder named Sources and
        clone the working_mesa_directory there.  You should also copy any other things from your existing .bash_profile file into
        this file so that you preserve your normal environment.  Make this your default .bash_profile file.
    3.  Under normal conditions, the .bash_profile file will get executed as soon as you enter the bash shell.  The first time you set
        it up you either need to source it or crank up a new terminal after making the changes in 2.  When .bash_profile executes, it
        will call other scripts which will ask you questions.  These questions need to be answered or the default values will be
        used.  There are three scripts are scripts are, general environment.sh, compilers.sh and mesa.sh.   These are located in the
        scripts folder.  There are also scripts in that folder for ifort, g95, gfortran, g77, f77 and pgi.  Things have been tested
        extensively using the Intel and PGI compilers.  It appears Gfortran also works, but I have not dome as much testing as
        withe the other two compilers.

        Here are the questions.  They are a bit different on the mac.  This is the description for a linux machine.
   
           * First question is location of sources - if you hit return then you get the default.
             > A line is printed giving you the OS and a Machine parameter
           * Second question is whether you are using MiC's - Thats the Intel Phi coprocessor. The default is not to use the MiC's
             If you hit return it defaults to none. 
           * Third question is what compiler are you using and whether you want 32 or 64 bit integers.  You need to explicitly say
             ifort and i8 if you want 64 bit integers.  The default is ifort and 32 bit integers.
           * Fourth question is about preprocessor flags.  In most cases take the default.
           * Fifth question is about MPI version.  You can use either intelmpi or mvapich2.
           * Sixth question is about linking model.  It can be static or dynamic.  Either is fine but you need to specify.
           * Seventh question is about layer option.  The answer depends on whether you are using threads or not.
           * Eighth question is about the Intel home directory.  This is where the ifort compiler and mkl libraries are located.  
             There is a default but that is something you need to set by hand in the script.
           * Ninth question asks if you want to clone mesa to your machine.  Normally you would have it but this would allow you to clone
             it from bitbucket.

                  a. There are defaults set for these questions
                  b. Best approach is to run the .bash_profile, answer the questions and look at the environment variables.

           * echo $MESA_HOME, $MESA_BIN, $MESA_TMP, $F90, $MESA_RUN, F90FLAGS and $MD_LIB to see some of the essential flags.
        
              I.   $MESA_HOME =  Location of your copy of mesa sources  
              II.  $MESA_BIN  =  Location of mesa executables  
              III. $MESA_TMP  = The directory for the wf, rint, int, etc. files.   
                   It should point to a place where you have plenty of room such as a large disk or SSD. 
                   It is essential to have a directory of this sort for each user, and perhaps even for 
                   each user on each machine that they use so that one users rwf file wont accidentally 
                   overwrite anothers.  While overwriting file names can be prevented by setting unique 
                   names, the default name is fixed.  Note also, wWhen you set this, dont forget to 
                   actually create the /xxx/user/tmp ddirectory.  I usually have a run directory 
                   on each machine and create a soft link to a directory with lots of space on that machine.  
                   One certainly does not want MESA_TMP to be physically on another machine as one can 
                   imagine the time it will take to move large files to and from that space.  
                   It has been shown useful for each user to set up a small directory in which 
                   they run mesa.  For example, I have a /xxx/mesa_working/run where I run the code and also create 
                   the softlink to the MESA_TMP directory.  There are subdirectories of this directory
                   such as  inp/ out/ chk/ to store needed input,output and check files.  
                   I also keep my runmesa script there which runs the mesa jobs.  If run in the 
                   background, there is a log file produced with timing statistics. 
              IV.  A $MESA_RUN directory should be set up where input files reside.  The run directory should be set up
                   on each machine.

## Makefiles.

    A lot of the original way in which the Makefiles were set up are now moot.
    In the main mesa directory there is a Makefile and a Makefile.inc  The Makefile is complex as it was 
    designed to do many tasks depending on what links and library files you want to make or clean.  
    Please look at it and become familiar with the various usages.  The Makefile.inc contains the  
    definitions of where libraries are located and commands needed to compile and load the links.  It is easy
    to compile and link now.  The makefiles for the individual links are much simplified and very few commands
    are needed to actually do things.
 
## Who do I talk to?

    * Barry Schneider
                * Email: bis@nist.gov
                * Phone: 301-975-4685