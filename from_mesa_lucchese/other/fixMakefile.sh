#
for SRC in m[0-9]* ; do
cp ~/Desktop/dir_mesa77/$SRC/Makefile $SRC
done
for KSRC in k* ; do
SRC=${KSRC/k/}
cp ~/Desktop/dir_mesa77/$SRC/Makefile $KSRC
sed -e "s/\/$SRC\//\/$KSRC\//g" < $KSRC/Makefile >stuff
mv stuff $KSRC/Makefile
done

for MF in m[0-9]* k* ; do sed -e "s/^FC =/# FC =/g" <$MF/Makefile >stuff ; mv stuff $MF/Makefile  ; done
egrep -h "FC =" k*/Makefile m[0-9]*/Makefile | sort -u

for MF in m[0-9]* k* ; do sed -e "s/^LD =/# LD =/g" <$MF/Makefile >stuff ; mv stuff $MF/Makefile  ; done
egrep -h "LD =" k*/Makefile m[0-9]*/Makefile | sort -u

for MF in m[0-9]* k* ; do sed -e "s/^LDFLAGS = /# LDFLAGS =/g" <$MF/Makefile >stuff ; mv stuff $MF/Makefile  ; done
for MF in m[0-9]* k* ; do sed -e "s/^LDFLAGS =/# LDFLAGS =/g" <$MF/Makefile >stuff ; mv stuff $MF/Makefile  ; done
egrep -h "LDFLAGS =" k*/Makefile m[0-9]*/Makefile | sort -u

for MF in m[0-9]* k* ; do sed -e "s/^FFLAGS = -c -g -finit-local-zero -g -finit-local-zero -fno-backslash/FFLAGS = \$(FFLAGSZB)/g"\
            <$MF/Makefile >stuff ; mv stuff $MF/Makefile  ; done
for MF in m[0-9]* k* ; do sed -e "s/^FFLAGS = -c -g -finit-local-zero -g -finit-local-zero/FFLAGS = \$(FFLAGSZ)/g"\
            <$MF/Makefile >stuff ; mv stuff $MF/Makefile  ; done
for MF in m[0-9]* k* ; do sed -e "s/^FFLAGS = -g -finit-local-zero -c/FFLAGS = \$(FFLAGSZ)/g"\
            <$MF/Makefile >stuff ; mv stuff $MF/Makefile  ; done
for MF in m[0-9]* k* ; do sed -e "s/^FFLAGS = -finit-local-zero -c -g /FFLAGS = \$(FFLAGSZ)/g"\
            <$MF/Makefile >stuff ; mv stuff $MF/Makefile  ; done
for MF in m[0-9]* k* ; do sed -e "s/^FFLAGS = -c -g -finit-local-zero -fno-automatic/FFLAGS = \$(FFLAGSZA)/g" \
            <$MF/Makefile >stuff ; mv stuff $MF/Makefile  ; done
for MF in m[0-9]* k* ; do sed -e "s/^FFLAGS = -c -g -finit-local-zero -O/FFLAGS = \$(FFLAGSOZ)/g" \
            <$MF/Makefile >stuff ; mv stuff $MF/Makefile  ; done
for MF in m[0-9]* k* ; do sed -e "s/^FFLAGS = -c -g -finit-local-zero /FFLAGS = \$(FFLAGSZ)/g" \
            <$MF/Makefile >stuff ; mv stuff $MF/Makefile  ; done
for MF in m[0-9]* k* ; do sed -e "s/^FFLAGS = -c -g -finit-local-zero/FFLAGS = \$(FFLAGSZ)/g"\
            <$MF/Makefile >stuff ; mv stuff $MF/Makefile  ; done
for MF in m[0-9]* k* ; do sed -e "s/^FFLAGS = -c -fno-backslash -finit-local-zero/FFLAGS = \$(FFLAGSZB)/g"\
            <$MF/Makefile >stuff ; mv stuff $MF/Makefile  ; done
for MF in m[0-9]* k* ; do sed -e "s/^FFLAGS = -c -finit-local-zero -g -finit-local-zero -O/FFLAGS = \$(FFLAGSOZ)/g"\
            <$MF/Makefile >stuff ; mv stuff $MF/Makefile  ; done
for MF in m[0-9]* k* ; do sed -e "s/^FFLAGS = -c -finit-local-zero -fno-backslash/FFLAGS = \$(FFLAGSZB)/g"\
            <$MF/Makefile >stuff ; mv stuff $MF/Makefile  ; done
for MF in m[0-9]* k* ; do sed -e "s/^FFLAGS = -c -finit-local-zero -finit-local-zero/FFLAGS = \$(FFLAGSZ)/g"\
            <$MF/Makefile >stuff ; mv stuff $MF/Makefile  ; done
for MF in m[0-9]* k* ; do sed -e "s/^FFLAGS = -c -finit-local-zero /FFLAGS = \$(FFLAGSZ)/g"\
            <$MF/Makefile >stuff ; mv stuff $MF/Makefile  ; done
for MF in m[0-9]* k* ; do sed -e "s/^FFLAGS = -c -O5 -g -finit-local-zero /FFLAGS = \$(FFLAGSO5Z)/g"\
            <$MF/Makefile >stuff ; mv stuff $MF/Makefile  ; done
for MF in m[0-9]* k* ; do sed -e "s/^FFLAGS = -c -O -finit-local-zero -g -finit-local-zero/FFLAGS = \$(FFLAGSOZ)/g"\
            <$MF/Makefile >stuff ; mv stuff $MF/Makefile  ; done
for MF in m[0-9]* k* ; do sed -e "s/^FFLAGS = -c -O -finit-local-zero -finit-local-zero/FFLAGS = \$(FFLAGSOZ)/g"\
            <$MF/Makefile >stuff ; mv stuff $MF/Makefile  ; done
for MF in m[0-9]* k* ; do sed -e "s/^FFLAGS = -c -O -finit-local-zero /FFLAGS = \$(FFLAGSOZ)/g"\
            <$MF/Makefile >stuff ; mv stuff $MF/Makefile  ; done
for MF in m[0-9]* k* ; do sed -e "s/^FFLAGS = -c -O /FFLAGS = \$(FFLAGSO)/g"\
            <$MF/Makefile >stuff ; mv stuff $MF/Makefile  ; done
for MF in m[0-9]* k* ; do sed -e "s/^FFLAGS = -c -O/FFLAGS = \$(FFLAGSO)/g"\
            <$MF/Makefile >stuff ; mv stuff $MF/Makefile  ; done

egrep -h "FFLAGS =" k*/Makefile m[0-9]*/Makefile | sort -u

for MF in m[0-9]* k* ; do sed -e "s/\/System\/Library\/Frameworks\/vecLib.framework\/versions\/A\/vecLib/\$(BLASUSE)/g"\
            <$MF/Makefile >stuff ; mv stuff $MF/Makefile  ; done

for MF in m[0-9]* k* ; do sed -e "s/\/System\/Library\/Frameworks\/Accelerate.framework\/Versions\/A\/Accelerate/\$(BLASUSEB)/g"\
            <$MF/Makefile >stuff ; mv stuff $MF/Makefile  ; done

egrep -h "BLAS *=|BLASLIB *=" k*/Makefile m[0-9]*/Makefile | sort -u
