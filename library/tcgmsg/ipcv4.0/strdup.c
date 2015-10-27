/* $Header: /msrc/proj/mss/tcgmsg/ipcv4.0/strdup.c,v 1.1 1994/02/23 16:55:19 d3g681 Exp $ */

#if defined(ULTRIX) || defined(SGI) || defined(NEXT) || defined(HPUX) \
                    || defined(DECOSF)
extern void *malloc();
#else
extern char *malloc();
#endif
extern char *strcpy();

char *strdup(s)
    char *s;
{
  char *new;

  if (new = malloc((unsigned) (strlen(s)+1)))
    (void) strcpy(new,s);

  return new;
}
