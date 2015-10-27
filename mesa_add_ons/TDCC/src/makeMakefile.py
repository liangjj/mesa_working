#!/usr/bin/env python

import os
import sys

makefilehead = open('Makefile.head').readlines()

sys.stdout.writelines(makefilehead)

print '### dependencies as found by makeMakefile.py'

f90files = {}
for fi in os.listdir('.'):
    fispl = os.path.splitext(fi)
    ext = fispl[1]
    if (ext == '.f90' or ext == '.F90' or ext == '.f95' or ext == '.F95'):
        f90files[fispl[0]] = ext

f90roots = f90files.keys()
f90roots.sort()

f90modules = {}
f90programs = {}

f90depends = {}

for firoot in f90roots:
    fi = firoot + f90files[firoot]
    f90depends[firoot] = []
    for li in open(fi,'r').readlines():
        li = li.lower()

        # remove comments from source
        pos = li.find('!')
        if (pos != -1):
            li = li[:pos]

        # remove strings from source - might contain "use" and confuse program
        try:
            while 1:
                pos = li.index("'")
                endpos = pos + 1 + li[pos+1:].find("'")
                li = li[:pos]+li[endpos:]

                pos = li.index('"')                
                endpos = pos + 1 + li[pos+1:].find("'")
                li = li[:pos]+li[endpos+1:]
        except ValueError:
            pass

        # check if this file contains a program
        try:
            li.index('program')
            f90programs[firoot] = 1
        except ValueError:
            pass

        # check if this file contains a module
        try:
            li.index('module')
            f90modules[firoot]  = 1
        except ValueError:
            pass

        # find modules used in use statements
        lil = li.split()
        try:
            while 1:
                pos = lil.index('use')
                module = lil[pos+1]
                if module[-1] == ';': module = module[:-1]
                if module[-1] == ',': module = module[:-1]
                f90depends[firoot].append(module)
                lil = lil[pos+2:]
        except ValueError:
            pass
    

for firoot in f90files.keys():
    linkobjs = {}
    tmpdic = {}
    for mo in f90depends[firoot]:
        if mo in f90roots:
            tmpdic[mo] = 1
            linkobjs[mo] = 1
        else:
            sys.stderr.write('did not find source file for module ' + mo + ', assuming it is somewhere else\n')        

    f90depends[firoot] = tmpdic.keys()
    f90depends[firoot].sort()

    linkobjs = linkobjs.keys()
    # also add all modules to linking that are used by modules used in this package
    linkschange = 1
    while linkschange:
        linkschange = 0
        for mo in linkobjs:
            # if there are dependencies for that link object
            if mo in f90depends.keys():
                for mo2 in f90depends[mo]:
                    # add file if we have the corresponding source file and it's not already in our list
                    if mo2 in f90roots and mo2 not in linkobjs:
                        linkobjs.append(mo2)
                        linkschange = 1
    linkobjs.sort()

    #try:
    #    linkobjs.index('tdcc_it')
    #    linkobjs.append('')
    #except ValueError:
    #    pass
    
    #if firoot == 'TDCC_IT':
    #    linkobjs.append('math_subroutines')

    if (len(f90depends[firoot]) > 0):
        objfi = firoot+'.o'
        if firoot in f90modules.keys():
            objfi += ' '+firoot+'.mod'
        print objfi, ':',
        for mo in f90depends[firoot]:
            print mo+'.mod',
        print
        print
        if firoot in f90programs.keys():
            exefi = firoot+'.exe'
            print exefi, ':',
            for mo in linkobjs:
                print mo+'.o',
            print
            print

