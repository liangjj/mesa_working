#
# the library makefiles need to initiall be modified by hand to remove dependence on the mesalib.a and to force
# creation
liblist="chr clams dpintgrl io math mdutil/sun sym util intgrl"

for MF in $liblist ; do sed -e "s/^FC =/# FC =/g" <$MF/Makefile >stuff ; mv stuff $MF/Makefile  ; done
for MF in $liblist ; do egrep -h "FC =|FC=" $MF/Makefile ; done | sort -u

for MF in $liblist ; do sed -e "s/^FFLAGS = -finit-local-zero -g -c -O2/FFLAGS = \$(FFLAGSO2Z)/g"\
            <$MF/Makefile >stuff ; mv stuff $MF/Makefile  ; done

for MF in $liblist ; do sed -e "s/^FFLAGS = -c -O2/FFLAGS = \$(FFLAGSO2)/g"\
            <$MF/Makefile >stuff ; mv stuff $MF/Makefile  ; done

for MF in $liblist ; do sed -e "s/^FFLAGS = -c/FFLAGS = \$(FFLAGSC)/g"\
            <$MF/Makefile >stuff ; mv stuff $MF/Makefile  ; done

for MF in $liblist ; do egrep -h "FFLAGS =" $MF/Makefile ; done | sort -u


