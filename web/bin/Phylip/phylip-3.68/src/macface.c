
#include "macface.h"

/* version 3.6 Copyright 1997-2004 by the University of Washington.
 * mac interfacing
 * this is only included and used for macs/powermacs
 * and needs to be used in conjunction with the customized LAMARC-library
 * which is a modification of the Metrowerks custom library.
 * fixmacfile() was written by Sean Lamont 1994.
 *
 * Peter Beerli 1997
 * beerli@genetics.washington.edu
 */

#include <SIOUX.h>
#include "SIOUXGlobals.h"
#include "SIOUXMenus.h"
#include "SIOUXWindows.h"

#include <console.h>
#include <signal.h>
#include <stdlib.h>
#include <stdio.h>
#include <Fonts.h> /*//bkoz*/
#include <Dialogs.h>
#include <AppleEvents.h>
#include <ToolUtils.h>
#ifndef OSX_CARBON
#include <DiskInit.h>
#include <Traps.h>
#endif
/*//defined in unix.c ...*/
#ifdef __cplusplus
extern "C" {
#endif

int  __system7present(void);
long __myraise(long signal);


#ifdef __cplusplus
}
#endif
#include <Files.h>
#include "phylip.h"

void fixmacfile(char *filename)
{
  FInfo  fndrinfo;
#ifdef OSX_CARBON
  FSSpec fileSpec;
  FSRef fileRef;

  FSPathMakeRef((unsigned char *)filename, &fileRef, NULL);
  FSGetCatalogInfo(&fileRef, kFSCatInfoNone, NULL, NULL, &fileSpec, NULL);
  FSpGetFInfo(&fileSpec, &fndrinfo);
  fndrinfo.fdType='TEXT';
  fndrinfo.fdCreator='MSWD';
  FSpSetFInfo(&fileSpec, &fndrinfo);

#else
  char filename1[FNMLNGTH];
  char filename2[FNMLNGTH];
  strcpy(filename1,filename);
  strcpy(filename2,filename);
  getfinfo(filename1,0,&fndrinfo);
  fndrinfo.fdType='TEXT';
  fndrinfo.fdCreator='MSWD';
  setfinfo(filename2,0,&fndrinfo);
#endif
}
/* new fixmacfile from interface.c is more versatile and may be desireable here
LM 7/27*/
/*
void fixmacfile(char *filename, long type, long creator)
{
FInfo  fndrinfo;
char filename1[FNMLNGTH];
char filename2[FNMLNGTH];
strcpy(filename1,filename);
strcpy(filename2,filename);
getfinfo(filename1,0,&fndrinfo);
fndrinfo.fdType=type;
fndrinfo.fdCreator=creator;
setfinfo(filename2,0,&fndrinfo);
}
*/



