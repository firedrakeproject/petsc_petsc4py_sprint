#ifndef lint
static char vcid[] = "$Id: fwd.c,v 1.8 1997/02/22 02:23:29 bsmith Exp balay $";
#endif
/*
      Code for manipulating files.
*/
#include "src/sys/src/files.h"

#undef __FUNC__  
#define __FUNC__ "PetscGetWorkingDirectory" /* ADIC Ignore */
/*@C
   PetscGetWorkingDirectory - Gets the current working directory.

   Input Paramters:
.  path - use to hold the result value
.  len  - maximum length of path

.keywords, system, get, current, working, directory
@*/
int PetscGetWorkingDirectory( char *path,int len )
{
#if defined(PARCH_sun4)
  getwd( path );
#elif defined(PARCH_nt)
  _getcwd( path, len );
#else
  getcwd( path, len );
#endif
  return 0;
}
