

int ftocstr(ds, dl, ss, sl)
    char *ds, *ss;      /* dst, src ptrs */
    int dl;             /* dst max len */
    int sl;             /* src len */
{
    char *p;

    for (p = ss + sl; --p >= ss && *p == ' '; ) ;
    sl = p - ss + 1;
    dl--;
    ds[0] = 0;
    if (sl > dl)
        return 1;
    while (sl--)
	(*ds++ = *ss++);
    *ds = '\0';
    return 0;
}


int ctofstr(ds, dl, ss)
	char *ds;		/* dest space */
	int dl;			/* max dest length */
	char *ss;		/* src string (0-term) */
{
    while (dl && *ss) {
	*ds++ = *ss++;
	dl--;
    }
    while (dl--)
	*ds++ = ' ';
    return 0;
}
