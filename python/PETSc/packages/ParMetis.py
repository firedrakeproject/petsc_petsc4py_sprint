from __future__ import generators
import config.base

class Configure(config.base.Configure):
  def __init__(self, framework):
    config.base.Configure.__init__(self, framework)
    self.headerPrefix = ''
    self.substPrefix  = ''
    self.found        = 0
    self.compilers    = self.framework.require('config.compilers', self)
    self.libraries    = self.framework.require('config.libraries', self)
    return

  def __str__(self):
    if self.found:
      desc = ['ParMetis:']	
      desc.append('  Version: '+self.version)
      desc.append('  Includes: '+str(self.include))
      desc.append('  Library: '+str(self.lib))
      return '\n'.join(desc)+'\n'
    else:
      return ''

  def setupHelp(self, help):
    import nargs
    help.addArgument('ParMetis', '-with-parmetis=<bool>',                nargs.ArgBool(None, 0, 'Activate ParMetis'))
    help.addArgument('ParMetis', '-with-parmetis-dir=<root dir>',        nargs.ArgDir(None, None, 'Specify the root directory of the ParMetis installation'))
    help.addArgument('ParMetis', '-with-parmetis-include=<dir>',         nargs.ArgDir(None, None, 'The directory containing metis.h'))
    help.addArgument('ParMetis', '-with-parmetis-lib=<lib>',             nargs.Arg(None, None, 'The ParMetis library or list of libraries'))
    help.addArgument('ParMetis', '-download-parmetis=<no,yes,ifneeded>', nargs.ArgFuzzyBool(None, 2, 'Install MPICH to provide ParMetis'))
    return

  def checkLib(self, libraries):
    '''Check for ParMetis_Init in libraries, which can be a list of libraries or a single library'''
    if not isinstance(libraries, list): libraries = [libraries]
    oldLibs = self.framework.argDB['LIBS']
    found   = self.libraries.check(libraries, 'ParMetis_Init')
    self.framework.argDB['LIBS'] = oldLibs
    return found

  def checkInclude(self, includeDir):
    '''Check that parmetis.h is present'''
    oldFlags = self.framework.argDB['CPPFLAGS']
    for inc in includeDir:
      self.framework.argDB['CPPFLAGS'] += ' -I'+inc
    found = self.checkPreprocess('#include <parmetis.h>\n')
    self.framework.argDB['CPPFLAGS'] = oldFlags
    return found

  def includeGuesses(self, path):
    '''Return all include directories present in path or its ancestors'''
    while path:
      dir = os.path.join(path, 'include')
      if os.path.isdir(dir):
        yield [dir]
      path = os.path.dirname(path)
    return

  def libraryGuesses(self, root = None):
    '''Return standard library name guesses for a given installation root'''
    if root:
      yield [os.path.join(root, 'lib', 'libparmetis.a'), os.path.join(root, 'lib', 'libmetis.a')]
    else:
      yield ['']
      yield ['parmetis','metis']
    return

  def generateGuesses(self):
    if 'with-parmetis-lib' in self.framework.argDB and 'with-parmetis-dir' in self.framework.argDB:
      raise RuntimeError('You cannot give BOTH ParMetis library with --with-parmetis-lib=<lib> and search directory with --with-parmetis-dir=<dir>')
    if self.framework.argDB['download-parmetis'] == 1:
      (name, lib, include) = self.downloadParMetis()
      yield (name, lib, include)
      raise RuntimeError('Downloaded ParMetis could not be used. Please check install in '+os.path.dirname(include[0])+'\n')
    # Try specified library and include
    if 'with-parmetis-lib' in self.framework.argDB:
      libs = self.framework.argDB['with-parmetis-lib']
      if not isinstance(libs, list): libs = [libs]
      if 'with-parmetis-include' in self.framework.argDB:
        includes = [[self.framework.argDB['with-parmetis-include']]]
      else:
        includes = self.includeGuesses('\n'.join(map(lambda inc: os.path.dirname(os.path.dirname(inc)), libs)))
      yield ('User specified library and includes', [libs], includes)
      raise RuntimeError('You set a value for --with-parmetis-lib, but '+str(self.framework.argDB['with-parmetis-lib'])+' cannot be used.\n')
    # Try specified installation root
    if 'with-parmetis-dir' in self.framework.argDB:
      dir = self.framework.argDB['with-parmetis-dir']
      if not (len(dir) > 2 and dir[1] == ':'):
        dir = os.path.abspath(dir)
      yield ('User specified installation root', self.libraryGuesses(dir), [[os.path.join(dir, 'include')]])
      raise RuntimeError('You set a value for --with-parmetis-dir, but '+self.framework.argDB['with-parmetis-dir']+' cannot be used.\n')
    # May not need to list anything
    yield ('Default compiler locations', self.libraryGuesses(), [[]])
    # Try configure package directories
    dirExp = re.compile(r'((p|P)ar)?(m|M)etis(-.*)?')
    for packageDir in self.framework.argDB['package-dirs']:
      packageDir = os.path.abspath(packageDir)
      if not os.path.isdir(packageDir):
        raise RuntimeError('Invalid package directory: '+packageDir)
      for f in os.listdir(packageDir):
        dir = os.path.join(packageDir, f)
        if not os.path.isdir(dir):
          continue
        if not dirExp.match(f):
          continue
        yield ('Package directory installation root', self.libraryGuesses(dir), [[os.path.join(dir, 'include')]])
    # Try /usr/local
    dir = os.path.abspath(os.path.join('/usr', 'local'))
    yield ('Frequent user install location (/usr/local)', self.libraryGuesses(dir), [[os.path.join(dir, 'include')]])
    # Try /usr/local/*parmetis*
    ls = os.listdir(os.path.join('/usr','local'))
    for dir in ls:
      if not dirExp.match(f):
        continue
      dir = os.path.join('/usr','local',dir)
      if os.path.isdir(dir):
        yield ('Frequent user install location (/usr/local/parmetis*)', self.libraryGuesses(dir), [[os.path.join(dir, 'include')]])
    # Try ~/parmetis*
    ls = os.listdir(os.getenv('HOME'))
    for dir in ls:
      if not dirExp.match(f):
        continue
      dir = os.path.join(os.getenv('HOME'),dir)
      if os.path.isdir(dir):
        yield ('Frequent user install location (~/parmetis*)', self.libraryGuesses(dir), [[os.path.join(dir, 'include')]])
    # Try PETSc location
    PETSC_DIR  = None
    PETSC_ARCH = None
    if 'PETSC_DIR' in self.framework.argDB and 'PETSC_ARCH' in self.framework.argDB:
      PETSC_DIR  = self.framework.argDB['PETSC_DIR']
      PETSC_ARCH = self.framework.argDB['PETSC_ARCH']
    elif os.getenv('PETSC_DIR') and os.getenv('PETSC_ARCH'):
      PETSC_DIR  = os.getenv('PETSC_DIR')
      PETSC_ARCH = os.getenv('PETSC_ARCH')
    if PETSC_ARCH and PETSC_DIR:
      pass
##      try:
##        libArgs = config.base.Configure.executeShellCommand('cd '+PETSC_DIR+'; make BOPT=g_c++ getmpilinklibs', log = self.framework.log)[0].strip()
##        incArgs = config.base.Configure.executeShellCommand('cd '+PETSC_DIR+'; make BOPT=g_c++ getmpiincludedirs', log = self.framework.log)[0].strip()
##        libArgs = self.splitLibs(libArgs)
##        incArgs = self.splitIncludes(incArgs)
##        if libArgs and incArgs:
##          yield ('PETSc location', [libArgs], [incArgs])
##      except RuntimeError:
##        # This happens with older Petsc versions which are missing those targets
##        pass
    # If necessary, download ParMetis
    if not self.found and self.framework.argDB['download-parmetis'] == 2:
      (name, lib, include) = self.downloadParMetis()
      yield (name, lib, include)
      raise RuntimeError('Downloaded ParMetis could not be used. Please check in install in '+os.path.dirname(include)+'\n')
    return

  def getDir(self):
    '''Find the directory containing ParMetis'''
    packages  = os.path.join(self.framework.argDB['PETSC_DIR'], 'packages')
    if not os.path.isdir(packages):
      os.mkdir(packages)
      self.framework.actions.addArgument('PETSc', 'Directory creation', 'Created the packages directory: '+packages)
    parmetisDir = None
    for dir in os.listdir(packages):
      if dir.startswith('ParMetis') and os.path.isdir(os.path.join(packages, dir)):
        parmetisDir = dir
    if parmetisDir is None:
      self.framework.logPrint('Could not locate already downloaded ParMetis')
      raise RuntimeError('Error locating ParMetis directory')
    return os.path.join(packages, mpichDir)

  def downloadParMetis(self):
    self.framework.logPrint('Downloading ParMetis')
    try:
      mpichDir = self.getDir()
      self.framework.logPrint('ParMetis already downloaded, no need to ftp')
    except RuntimeError:
      import urllib

      packages = os.path.join(self.framework.argDB['PETSC_DIR'], 'packages')
      try:
        urllib.urlretrieve('ftp://ftp.mcs.anl.gov/pub/petsc/parmetis.tar.gz', os.path.join(packages, 'parmetis.tar.gz'))
      except Exception, e:
        raise RuntimeError('Error downloading ParMetis: '+str(e))
      try:
        config.base.Configure.executeShellCommand('cd packages; gunzip parmetis.tar.gz', log = self.framework.log)
      except RuntimeError, e:
        raise RuntimeError('Error unzipping parmetis.tar.gz: '+str(e))
      try:
        config.base.Configure.executeShellCommand('cd packages; tar -xf parmetis.tar', log = self.framework.log)
      except RuntimeError, e:
        raise RuntimeError('Error doing tar -xf parmetis.tar: '+str(e))
      os.unlink(os.path.join(packages, 'parmetis.tar'))
      self.framework.actions.addArgument('ParMetis', 'Download', 'Downloaded ParMetis into '+self.getDir())
    # Get the ParMetis directories
    parmetisDir = self.getDir()
    installDir = os.path.join(parmetisDir, self.framework.argDB['PETSC_ARCH'])
    if not os.path.isdir(installDir):
      os.mkdir(installDir)
    # Configure and Build ParMetis
    args = ['--prefix='+installDir, '-cc='+self.framework.argDB['CC']]
    args = ' '.join(args)
    try:
      fd      = file(os.path.join(installDir,'config.args'))
      oldargs = fd.readline()
      fd.close()
    except:
      oldargs = ''
    if not oldargs == args:
      self.framework.log.write('Have to rebuild ParMetis oldargs = '+oldargs+' new args '+args+'\n')
      try:
        output  = config.base.Configure.executeShellCommand('cd '+parmetisDir+';./configure '+args, timeout=900, log = self.framework.log)[0]
      except RuntimeError, e:
        raise RuntimeError('Error running configure on ParMetis: '+str(e))
      try:
        output  = config.base.Configure.executeShellCommand('cd '+parmetisDir+';make; make install', timeout=2500, log = self.framework.log)[0]
      except RuntimeError, e:
        raise RuntimeError('Error running make; make install on ParMetis: '+str(e))
      fd = file(os.path.join(installDir,'config.args'), 'w')
      fd.write(args)
      fd.close()
      self.framework.actions.addArgument('MPI', 'Install', 'Installed ParMetis into '+installDir)
    lib     = [[os.path.join(installDir, 'lib', 'libparmetis.a'), os.path.join(installDir, 'lib', 'libmetis.a')]]
    include = [[os.path.join(installDir, 'include')]]
    return ('Downloaded ParMetis', lib, include)

  def configureVersion(self):
    '''Determine the ParMetos version'''
    raise RuntimeError('I have no fucking clue how to get the version')
    if self.framework.argDB['can-execute']:
      output, status = self.outputMPIRun('#include <stdio.h>\n#include <mpi.h>\n', 'int ver, subver;\n if (MPI_Get_version(&ver, &subver));\nprintf("%d.%d\\n", ver, subver)\n')
      if not status:
        # need to strip out information from batch system
        f = re.match('([0-9]*.[0-9]*)',output)
        if not f: return 'Unknown'
        return f.group()
    return 'Unknown'

  def configureLibrary(self):
    '''Find all working MPI installations and then choose one'''
    functionalParMetis = []
    for (name, libraryGuesses, includeGuesses) in self.generateGuesses():
      self.framework.logPrint('================================================================================')
      self.framework.logPrint('Checking for a functional ParMetis in '+name)
      self.lib     = None
      self.include = None
      found        = 0
      for libraries in libraryGuesses:
        if self.checkLib(libraries):
          self.lib = libraries
          for includeDir in includeGuesses:
            if self.checkInclude(includeDir):
              self.include = includeDir
              found = 1
              break
          if found:
            break
      if not found: continue
      version = self.executeTest(self.configureVersion)
      self.found = 1
      functionalParMetis.append((name, self.lib, self.include, version))
      if not self.framework.argDB['with-alternatives']:
        break
    # User chooses one or take first (sort by version)
    if self.found:
      self.name, self.lib, self.include, self.version = functionalParMetis[0]
      self.framework.logPrint('Choose MPI '+self.version+' in '+self.name)
    else:
      self.framework.logPrint('Could not locate any functional ParMetis')
    return

  def setOutput(self):
    '''Add defines and substitutions
       - HAVE_PARMETIS is defined if a working ParMetis is found
       - PARMETIS_INCLUDE and PARMETIS are command line arguments for the compile and link
       - PARMETIS_INCLUDE_DIR is the directory containing metis.h
       - PARMETIS_LIBRARY is the list of ParMetis libraries'''
    if self.found:
      self.addDefine('HAVE_PARMETIS', 1)
      if self.include:
        self.addSubstitution('PARMETIS_INCLUDE',     ' '.join(['-I'+inc for inc in self.include]))
        self.addSubstitution('PARMETIS_INCLUDE_DIR', self.include[0])
      else:
        self.addSubstitution('PARMETIS_INCLUDE',     '')
        self.addSubstitution('PARMETIS_INCLUDE_DIR', '')
      if self.lib:
        self.addSubstitution('PARMETIS_LIB',     ' '.join(map(self.libraries.getLibArgument, self.lib)))
        self.addSubstitution('PARMETIS_LIBRARY', self.lib)
      else:
        self.addSubstitution('PARMETIS_LIB',     '')
        self.addSubstitution('PARMETIS_LIBRARY', '')
    else:
      self.addSubstitution('PARMETIS_INCLUDE', '')
      self.addSubstitution('PARMETIS_LIB', '')
#    self.framework.packages.append(self)
    return

  def configure(self):
    if self.framework.argDB['with-parmetis']  and not self.framework.argDB['with-64-bit-ints']:
      self.executeTest(self.configureLibrary)
    self.setOutput()
    return
