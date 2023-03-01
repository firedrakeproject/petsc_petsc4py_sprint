#!/usr/bin/env python
#!/bin/env python
#
#    Generates fortran stubs for PETSc using the Sowing bfort program
#
#    This tool looks for the values MANSEC and [BFORT]SUBMANSEC (where BFORTSUBMANSEC has priority over SUBMANSEC)
#    defined in the makefile (or include file when there is no makefile)
#
#    The F90 generated interface fils are stored in src/MANSEC/f90-mod/ftn-auto-interfaces/petsc[BFORT]SUBMANSEC.h90
#
#    These are then included by the petsc[[BFORT]SUB]MANSECmod.F90 files to create the Fortran module files
#
#    SUBMANSEC (but not BFORTSUBMANSEC) is also used (in the documentation generating part of PETSc) to determine what
#    directory in doc/manualpages/ the manual pages are deposited.
#
#    An example of when the BFORTSUBMANSEC may be different than SUBMANSEC is src/dm/label/impls/ephemeral/plex where we would like
#    the documentation to be stored under DMLabel but the function interfaces need to go into the DMPLEX Fortran module
#    (not the DM Fortran module) since they depend on DMPlexTransform.
#
from __future__ import print_function
import os

import subprocess
try:
  from subprocess import check_output
except ImportError:
  def check_output(*popenargs, **kwargs):
    """Implementation from Python-2.7 subprocess.check_output for use with
    Python-2.6 which does not provide check_output.
    """
    if 'stdout' in kwargs:
      raise ValueError('stdout argument not allowed, it will be overridden.')
    process = subprocess.Popen(stdout=subprocess.PIPE, *popenargs, **kwargs)
    output, unused_err = process.communicate()
    retcode = process.poll()
    if retcode:
      cmd = kwargs.get("args")
      if cmd is None:
        cmd = popenargs[0]
      raise subprocess.CalledProcessError(retcode, cmd, output=output)
    return output

#
def FixFile(filename):
  ''' Fixes the C fortran stub files generated by bfort'''
  import re

  def findLineCol(filename, string):
    with open(filename) as f:
      for l, line in enumerate(f, 1):
        c = line.find(string) + 1
        if c > 0: return l, c
    return 0, 0

  # https://gitlab.com/petsc/petsc/-/issues/360#note_231598172
  l, c = findLineCol(filename, '\00')
  if l > 0: print('WARNING: Found null character in generated Fortran stub file:\n  %s:%d:%d\n' % (filename, l, c))

  with open(filename) as ff:
    data = ff.read()

  data = re.subn('\00','',data)[0]
  data = re.subn('\nvoid ','\nPETSC_EXTERN void ',data)[0]
  data = re.subn('\nPetscErrorCode ','\nPETSC_EXTERN void ',data)[0]
  data = re.subn('Petsc([ToRm]*)Pointer\(int\)','Petsc\\1Pointer(void*)',data)[0]
  data = re.subn('PetscToPointer\(a\) \(a\)','PetscToPointer(a) (*(PetscFortranAddr *)(a))',data)[0]
  data = re.subn('PetscFromPointer\(a\) \(int\)\(a\)','PetscFromPointer(a) (PetscFortranAddr)(a)',data)[0]
  data = re.subn('PetscToPointer\( \*\(int\*\)','PetscToPointer(',data)[0]
  data = re.subn('MPI_Comm comm','MPI_Comm *comm',data)[0]
  data = re.subn('\(MPI_Comm\)PetscToPointer\( \(comm\) \)','MPI_Comm_f2c(*(MPI_Fint*)(comm))',data)[0]
  data = re.subn('\(PetscInt\* \)PetscToPointer','',data)[0]
  data = re.subn('\(Tao\* \)PetscToPointer','',data)[0]
  data = re.subn('\(TaoConvergedReason\* \)PetscToPointer','',data)[0]
  data = re.subn('\(TaoLineSearch\* \)PetscToPointer','',data)[0]
  data = re.subn('\(TaoLineSearchConvergedReason\* \)PetscToPointer','',data)[0]
  match = re.compile(r"""\b(PETSC|TAO)(_DLL|VEC_DLL|MAT_DLL|DM_DLL|KSP_DLL|SNES_DLL|TS_DLL|FORTRAN_DLL)(EXPORT)""")
  data = match.sub(r'',data)

  with open(filename, 'w') as ff:
    ff.write('#include "petscsys.h"\n#include "petscfix.h"\n#include "petsc/private/fortranimpl.h"\n'+data)

  # https://gitlab.com/petsc/petsc/-/issues/360#note_231598172
  l, c = findLineCol(filename, '\00')
  if l > 0: print('WARNING: Found null character in generated Fortran stub file after generatefortranstubs.py processing:\n  %s:%d:%d\n' % (filename, l, c))

def FindSource(filename):
  import os.path
  gendir, fname = os.path.split(filename)
  base, ext = os.path.splitext(fname)
  sdir, ftn_auto = os.path.split(gendir)
  if ftn_auto != 'ftn-auto': return None # Something is wrong, skip
  sfname = os.path.join(sdir, base[:-1] + ext)
  return sfname
  sourcefile = FindSource(filename)
  if sourcefile and os.path.isfile(sourcefile):
    import shutil
    shutil.copystat(sourcefile, filename)
  return

def FixDir(petscdir,dir,verbose):
  ''' Fixes a directory of files generated by bfort.
      + Fixes the C stub files to be compilable C
      + Generates a makefile
      + copies over Fortran interface files that are generated'''
  import re

  submansec = 'unknown'
  mansec = 'unknown'
  bfortsubmansec = 'unknown'
  cnames = []
  hnames = []
  parentdir = os.path.abspath(os.path.join(dir,'..'))
  for f in os.listdir(dir):
    ext = os.path.splitext(f)[1]
    if ext == '.c' or ext == '.cxx':
      FixFile(os.path.join(dir, f))
      cnames.append(f)
    elif ext == '.h90':
      hnames.append(f)
  if (cnames != [] or hnames != []):
    mfile=os.path.abspath(os.path.join(parentdir,'makefile'))
    try:
      fd=open(mfile,'r')
    except:
      print('Error! missing file:', mfile)
      return
    inbuf = fd.read()
    fd.close()
    cppflags = ""
    libbase = ""
    for line in inbuf.splitlines():
      if line.find('CPPFLAGS') >=0:
        cppflags = line
      if line.find('LIBBASE') >=0:
        libbase = line
      elif line.find('SUBMANSEC') >=0:
        submansec = line.split('=')[1].lower().strip()
      elif line.find('BFORTSUBMANSEC') >=0:
        bfortsubmansec = line.split('=')[1].lower().strip()
      elif line.find('MANSEC') >=0:
        submansec = line.split('=')[1].lower().strip()
      if line.find('MANSEC') >=0 and not line.find('SUBMANSEC') >=0:
        mansec = line.split('=')[1].lower().strip()

    if not bfortsubmansec == 'unknown':
      submansec = bfortsubmansec

    # now assemble the makefile
    outbuf  =  '\n'
    outbuf +=  "#requiresdefine   'PETSC_HAVE_FORTRAN'\n"
    outbuf +=  'ALL: lib\n'
    outbuf +=   cppflags + '\n'
    outbuf +=  'CFLAGS   =\n'
    outbuf +=  'FFLAGS   =\n'
    outbuf +=  'SOURCEC  = ' +' '.join(cnames)+ '\n'
    outbuf +=  'SOURCEF  =\n'
    outbuf +=  'SOURCEH  = ' +' '.join(hnames)+ '\n'
    outbuf +=  'DIRS     =\n'
    outbuf +=  libbase + '\n'
    outbuf +=  'include ${PETSC_DIR}/lib/petsc/conf/variables\n'
    outbuf +=  'include ${PETSC_DIR}/lib/petsc/conf/rules\n'
    outbuf +=  'include ${PETSC_DIR}/lib/petsc/conf/test\n'

    ff = open(os.path.join(dir, 'makefile'), 'w')
    ff.write(outbuf)
    ff.close()

  # if dir is empty - remove it
  if os.path.exists(dir) and os.path.isdir(dir) and os.listdir(dir) == []:
    os.rmdir(dir)

  # save Fortran interface file generated (it is merged with others in a post-processing step)
  for filename in [f for f in os.listdir(parentdir) if re.match(r'f90module[0-9]+.f90', f)]:
    modfile = os.path.join(parentdir, filename)
    if os.path.exists(modfile):
      if verbose: print('Generating F90 interface for '+modfile)
      fd = open(modfile)
      txt = fd.read()
      fd.close()
      if txt:
        if not os.path.isdir(os.path.join(petscdir,'src',mansec,'f90-mod','ftn-auto-interfaces')): os.mkdir(os.path.join(petscdir,'src',mansec,'f90-mod','ftn-auto-interfaces'))
        if not os.path.isdir(os.path.join(petscdir,'src',mansec,'f90-mod','ftn-auto-interfaces',submansec+'-tmpdir')): os.mkdir(os.path.join(petscdir,'src',mansec,'f90-mod','ftn-auto-interfaces',submansec+'-tmpdir'))
        fname =  os.path.join(petscdir,'src',mansec,'f90-mod','ftn-auto-interfaces',submansec+'-tmpdir',os.path.relpath(parentdir,petscdir).replace('/','_')+'.h90')
        with open(fname,'a') as fd:
          fd.write(txt)
      os.remove(modfile)

def PrepFtnDir(dir):
  ''' Generate a ftn-auto directory if needed'''
  import shutil
  if os.path.exists(dir) and not os.path.isdir(dir):
    raise RuntimeError('Error - specified path is not a dir: ' + dir)
  elif not os.path.exists(dir):
    os.mkdir(dir)
  else:
    files = os.listdir(dir)
    for file in files:
      if os.path.isdir(os.path.join(dir,file)): shutil.rmtree(os.path.join(dir,file))
      else: os.remove(os.path.join(dir,file))
  return

def processDir(petscdir, bfort, verbose, dirpath, dirnames, filenames):
  ''' Runs bfort on a directory and then fixes the files generated by bfort including moving generated F90 fortran interface files'''
  outdir = os.path.join(dirpath,'ftn-auto')
  newls = [l for l in filenames if os.path.splitext(l)[1] in ['.c','.h','.cxx','.cu']]
  if newls:
    PrepFtnDir(outdir)
    options = ['-dir '+outdir, '-mnative', '-ansi', '-nomsgs', '-noprofile', '-anyname', '-mapptr',
               '-mpi', '-shortargname', '-ferr', '-ptrprefix Petsc', '-ptr64 PETSC_USE_POINTER_CONVERSION',
               '-fcaps PETSC_HAVE_FORTRAN_CAPS', '-fuscore PETSC_HAVE_FORTRAN_UNDERSCORE',
               '-f90mod_skip_header','-on_error_abort']
    split_ct = 10
    for i in range(0, len(newls), split_ct):
      cmd = 'BFORT_CONFIG_PATH='+os.path.join(petscdir,'lib','petsc','conf')+' '+bfort+' '+' '.join(options+newls[i:i+split_ct])+' -f90modfile f90module'+str(i)+'.f90'
      try:
        output = check_output(cmd, cwd=dirpath, shell=True, stderr=subprocess.STDOUT)
      except subprocess.CalledProcessError as e:
        raise SystemError(str(e)+'\nIn '+dirpath+'\n'+e.output.decode(encoding='UTF-8',errors='replace'));
    FixDir(petscdir,outdir,verbose)

  # remove from list of subdirectories all directories without source code
  dirnames[:] = [name for name in dirnames
                 if name not in ['SCCS', 'output', 'BitKeeper', 'examples', 'externalpackages', 'bilinear', 'ftn-auto','ftn-auto-interfaces', 'fortran','bin','maint','ftn-custom','config','f90-custom','ftn-kernels','linter']
                 and not name.startswith(".")
                 # skip for ./configure generated $PETSC_ARCH directories
                 and not os.path.isdir(os.path.join(dirpath,name,'lib','petsc'))
                 and not os.path.isdir(os.path.join(dirpath,name,'lib','petsc-conf'))
                 and not os.path.isdir(os.path.join(dirpath,name,'conf'))
                 # skip include/petsc directory
                 and name != 'petsc']
  return

def updatePetscTypesFromMansec(types, path):
  for file in os.listdir(path):
    if file.endswith('.h'):
      with open(os.path.join(path,file)) as fd:
        txtlst = fd.readlines()
        lsts = [l.strip().split(' ') for l in txtlst if ' type ' in l]
        # l[0] == ! means comment, don't include comments
        newTypes = set(l[l.index('type')+1] for l in lsts if '!' not in l[0])
        types.update(newTypes)
  return types

def checkHandWrittenF90Interfaces(badSrc, path):
  import re
  for file in os.listdir(path):
    if file.endswith('.h90') or file.endswith('.F90'):
      with open(os.path.join(path,file),'r') as fdr:
        lineno = 1
        raw = fdr.read()
        for ibuf in re.split('(?i)\n\s*interface',raw):
          res = re.search('(.*)(\s+end\s+interface)',ibuf,flags=re.DOTALL|re.IGNORECASE)
          try:
            lines = res.group(0).split('\n')
            useLine = [(s.strip(),idx+lineno,os.path.join(path,file)) for idx,s in enumerate(lines) if 'use petsc' in s and 'only:' not in s]
            badSrc.extend(useLine)
          except AttributeError:
            # when re.search comes up empty
            pass
          except IndexError:
            # "use" was not in res.group(0)
            pass
          lineno = lineno+ibuf.count('\n')
  return badSrc

def processf90interfaces(petscdir,verbose):
  import shutil
  ''' Takes all the individually generated fortran interface files and merges them into one for each mansec'''
  ptypes = set()
  mansecs = ['sys','vec','mat','dm','ksp','snes','ts','tao']
  mansecF90Dirs = [os.path.join(petscdir,'src',ms,'f90-mod') for ms in mansecs]
  badSrc = []
  for msfd in mansecF90Dirs:
    badSrc = checkHandWrittenF90Interfaces(badSrc, msfd)
    ptypes = updatePetscTypesFromMansec(ptypes,msfd)
  for src in badSrc:
    print('Importing entire package: "'+src[0]+'" line '+str(src[1])+' file '+src[2])
  if len(badSrc): raise RuntimeError
  for msfd in mansecF90Dirs:
    msfad = os.path.join(msfd,'ftn-auto-interfaces')
    for submansec in os.listdir(msfad):
      if verbose: print('Processing F90 interface for '+submansec)
      if os.path.isdir(os.path.join(msfad,submansec)):
        submansec = submansec[:-7]
        f90inc = os.path.join(msfad,'petsc'+submansec+'.h90')
        tmpDir = os.path.join(msfad,submansec+'-tmpdir')
        with open(f90inc,'w') as fd:
          for sfile in os.listdir(tmpDir):
            if verbose: print('  Copying in '+sfile)
            with open(os.path.join(tmpDir,sfile),'r') as fdr:
              for ibuf in fdr.read().split('      subroutine')[1:]:
                ibuf = '      subroutine'+ibuf
                ibuf = ibuf.replace('integer z','PetscErrorCode z')
                ibuf = ibuf.replace('integer a ! MPI_Comm','MPI_Comm a ! MPI_Comm')
                plist = [p for p in ptypes if ' '+p[1:]+' ' in ibuf]
                if plist: ibuf = ibuf.replace(')',')\n       import '+','.join(set(plist)),1)
                fd.write(ibuf)
        shutil.rmtree(tmpDir)
  # FixDir(petscdir,os.path.join(petscdir,'include','petsc','finclude','ftn-auto-interfaces'),verbose)
  return

def main(petscdir,bfort,dir,verbose):
  for dirpath, dirnames, filenames in os.walk(dir):
    fnames = [i for i in filenames if not i.find('#') > -1]
    processDir(petscdir, bfort, verbose, dirpath, dirnames, fnames)
  return
#
# generatefortranstubs bfortexectuable -verbose            -----  generates fortran stubs for a directory and all its children
# generatefortranstubs -merge  -verbose                    -----  merges fortran 90 interfaces definitions that have been generated
#
if __name__ ==  '__main__':
  import sys
  if len(sys.argv) < 2: sys.exit('Must give the BFORT program or -merge as the first argument')
  petscdir = os.environ['PETSC_DIR']
  if len(sys.argv) > 2: verbose = 1
  else: verbose = 0
  if not sys.argv[1].startswith('-'):
    main(petscdir,sys.argv[1],os.getcwd(),verbose)
  else:
    processf90interfaces(petscdir,verbose)
