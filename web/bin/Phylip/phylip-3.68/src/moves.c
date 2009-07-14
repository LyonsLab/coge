
#include "phylip.h"
#include "moves.h"


void inpnum(long *n, boolean *success)
{
  /* used by dnamove, dolmove, move, & retree */
  int fields;
  char line[100];

#ifdef WIN32
  phyFillScreenColor();
#endif
  fflush(stdout);
  getstryng(line);
  *n = atof(line);
  fields = sscanf(line,"%ld",n);
  *success = (fields == 1);
}  /* inpnum */


void prereverse(boolean ansi)
{
  /* turn on reverse video */
  printf(ansi ? "\033[7m": "");
}  /* prereverse */


void postreverse(boolean ansi)
{
  /* turn off reverse video */
  printf(ansi ? "\033[0m" : "");
}  /* postreverse */


void chwrite(Char ch, long num, long *pos, long leftedge, long screenwidth)
{
  long i;

  for (i = 1; i <= num; i++) {
    if ((*pos) >= leftedge && (*pos) - leftedge + 1 < screenwidth)
      putchar(ch);
    (*pos)++;
  }
}  /* chwrite */


void nnwrite(long nodenum,long num,long *pos,long leftedge,long screenwidth)
{
  long i, leftx;

  leftx = leftedge - (*pos);
  if ((*pos) >= leftedge && (*pos) - leftedge + num < screenwidth)
    printf("%*ld", (int)num, nodenum);
  else if (leftx > 0 && leftx < 3)
    for(i=0;i<num-leftx;i++)
      printf(" ");
  (*pos) += num;
}  /* nnwrite */


void stwrite(const char *s,long length,long *pos,long leftedge,
                long screenwidth)
{

  if ((*pos) >= leftedge && (*pos) - leftedge + 1 < screenwidth)
    printf("%*s", (int)length, s);
  (*pos) += length;
}  /* stwrite */


void help(const char *letters)
{
  /* display help information */
  char input[100];

  printf("\n\nR Rearrange a tree by moving a node or group\n");
  printf("# Show the states of the next %s that doesn't fit tree\n", letters);
  printf("+ Show the states of the next %s\n", letters);
  printf("-         ...     of the previous %s\n", letters);
  printf("S Show the states of a given %s\n", letters);
  printf(". redisplay the same tree again\n");
  printf("T Try all possible positions of a node or group\n");
  printf("U Undo the most recent rearrangement\n");
  printf("W Write tree to a file\n");
  printf("O select an Outgroup for the tree\n");
  printf("F Flip (rotate) branches at a node\n");
  printf("H Move viewing window to the left\n");
  printf("J Move viewing window downward\n");
  printf("K Move viewing window upward\n");
  printf("L Move viewing window to the right\n");
  printf("C show only one Clade (subtree) (useful if tree is too big)\n");
  printf("? Help (this screen)\n");
  printf("Q (Quit) Exit from program\n");
  printf("X Exit from program\n\n\n");
  printf("TO CONTINUE, PRESS ON THE Return OR Enter KEY");
#ifdef WIN32
  phyFillScreenColor();
#endif
  fflush(stdout);
  getstryng(input);
}  /* help */


void treeoptions(boolean waswritten, Char *ch, FILE **outtree,
                        Char *outtreename, Char *progname)
{ /* interactively get options for writing a tree */
  char input[100];

  if (waswritten) {
    printf("\nTree file already was open.\n");
    printf("   A   Add to this tree to tree file\n");
    printf("   R   Replace tree file contents by this tree\n");
    printf("   F   Write out tree to a different tree file\n");
    printf("   N   Do Not write out this tree\n");
    do {
      printf("Which should we do? ");
#ifdef WIN32
      phyFillScreenColor();
#endif
      fflush(stdout);
      getstryng(input);
      *ch  = input[0];
      uppercase(ch);
    } while (*ch != 'A' && *ch != 'R' && *ch != 'N' && *ch != 'F');
  }
  if (*ch == 'F'){
    outtreename[0] = '\0';
    while (outtreename[0] =='\0'){
      printf("Please enter a tree file name>");
#ifdef MAC
      fixmacfile(outtreename);
#endif
#ifdef WIN32
      phyFillScreenColor();
#endif
      fflush(stdout);
      getstryng(outtreename);
    }
    FClose(*outtree);
  }
  if (*ch == 'R' || *ch == 'A' || *ch == 'F' || !waswritten){
    openfile(outtree,outtreename,"output tree file",
                       (*ch == 'A' && waswritten) ? "a" : "w",
             progname,outtreename);
  }
}  /* treeoptions */


void window(adjwindow action, long *leftedge, long *topedge, long hscroll,
                        long vscroll, long treelines, long screenlines,
                        long screenwidth, long farthest, boolean subtree)
{
  /* move viewing window of tree */
  switch (action) {

  case left:
    if (*leftedge != 1)
      *leftedge -= hscroll;
    break;

  case downn:
    /* The 'topedge + 6' is needed to allow downward scrolling
       when part of the tree is above the screen and only 1 or 2 lines
       are below it. */
    if (treelines - *topedge + 6 >= screenlines)
      *topedge += vscroll;
    break;

  case upp:
    if (*topedge != 1)
      *topedge -= vscroll;
    break;

  case right:
    if ((farthest + 6 + nmlngth + ((subtree) ? 8 : 0)) >
        (*leftedge + screenwidth))
      *leftedge += hscroll;
    break;
  }
}  /* window */


void pregraph(boolean ansi)
{
  /* turn on graphic characters */
  /* used in move & dolmove */
  printf(ansi ? "\033(0" : "");
}  /* pregraph */


void pregraph2(boolean ansi)
{
  /* turn on graphic characters */
  /* used in dnamove & retree */
  if (ansi) {
    printf("\033(0");
    printf("\033[10m");
  }
}  /* pregraph2 */


void postgraph(boolean ansi)
{
  /* turn off graphic characters */
  /* used in move & dolmove */
  printf(ansi ? "\033(B" : "");
}  /* postgraph */


void postgraph2(boolean ansi)
{
  /* turn off graphic characters */
  /* used in dnamove & retree */
  if (ansi) {
    printf("\033[11m");
    printf("\033(B");
  }
}  /* postgraph2 */


void nextinc(long *dispchar, long *dispword, long *dispbit, long chars,
                        long bits, boolean *display, steptr numsteps, steptr weight)
{
  /* show next incompatible character */
  /* used in move & dolmove */
  long disp0;
  boolean done;

  *display = true;
  disp0 = *dispchar;
  done = false;
  do {
    (*dispchar)++;
    if (*dispchar > chars) {
      *dispchar = 1;
      done = (disp0 == 0);
    }
  } while (!(numsteps[*dispchar - 1] >
             weight[*dispchar - 1] ||
             *dispchar == disp0 || done));
  *dispword = (*dispchar - 1) / bits + 1;
  *dispbit = (*dispchar - 1) % bits + 1;
}  /* nextinc */


void nextchar(long *dispchar, long *dispword, long *dispbit, long chars,
                        long bits, boolean *display)
{
  /* show next character */
  /* used in move & dolmove */
  *display = true;
  (*dispchar)++;
  if (*dispchar > chars)
    *dispchar = 1;
  *dispword = (*dispchar - 1) / bits + 1;
  *dispbit = (*dispchar - 1) % bits + 1;
}  /* nextchar */


void prevchar(long *dispchar, long *dispword, long *dispbit, long chars,
                        long bits, boolean *display)
{
  /* show previous character */
  /* used in move & dolmove */
  *display = true;
  (*dispchar)--;
  if (*dispchar < 1)
    *dispchar = chars;
  *dispword = (*dispchar - 1) / bits + 1;
  *dispbit = (*dispchar - 1) % bits + 1;
}  /* prevchar */


void show(long *dispchar, long *dispword, long *dispbit, long chars,
                        long bits, boolean *display)
{
  /* used in move & dolmove */
  long i;
  boolean ok;

  do {
    printf("SHOW: (Character number or 0 to see none)? ");
    inpnum(&i, &ok);
    ok = (ok && (i == 0 || (i >= 1 && i <= chars)));
    if (ok && i != 0) {
      *display = true;
      *dispchar = i;
      *dispword = (i - 1) / bits + 1;
      *dispbit = (i - 1) % bits + 1;
    }
    if (ok && i == 0)
      *display = false;
  } while (!ok);
}  /* show */

