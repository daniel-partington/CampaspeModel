# Build the executables that are used in the flopy autotests
import os
import sys
import shutil
import platform
import subprocess
import flopy

try:
    import pymake
except:
    print('pymake is not installed...will not build executables')
    pymake = None

fc = 'gfortran'
cc = 'gcc'
dbleprec = False
# bindir should be in the user path to run flopy tests with appropriate
# executables
#
# by default bindir will be in user directory
# On windows will be C:\\Users\\username\\.local\\bin
# On linux and osx will be /Users/username/.local/bin
bindir = os.path.join(os.path.expanduser('~'), '.local', 'bin')
bindir = os.path.abspath(bindir)
# pass --bindir path/to/directory to define a different bin dir
for ipos, arg in enumerate(sys.argv):
    if arg.lower() == '--bindir':
        bindir = sys.argv[ipos + 1]
    elif arg.lower() == '--dryrun':
        print('will perform dryrun and not build executables')
        pymake = None
print(bindir)
if not os.path.exists(bindir):
    os.makedirs(bindir, exist_ok=True)

def test_build_mfnwt():
    if pymake is None:
        return
    starget = 'MODFLOW-NWT'
    exe_name = 'mfnwt'
    dirname = 'MODFLOW-NWT_1.1.3'
    url = "http://water.usgs.gov/ogw/modflow-nwt/{0}.zip".format(dirname)

    build_target(starget, exe_name, url, dirname)

    return

def run_cmdlist(cmdlist, cwd='.'):
    proc = subprocess.Popen(cmdlist, shell=False,
                            stdout=subprocess.PIPE,
                            stderr=subprocess.STDOUT,
                            cwd=cwd)
    stdout_data, stderr_data = proc.communicate()
    if proc.returncode != 0:
        if isinstance(stdout_data, bytes):
            stdout_data = stdout_data.decode('utf-8')
            stderr_data = stderr_datab.decode('utf-8')
        msg = '{} failed\n'.format(cmdlist) + \
              'status code:\n{}\n'.format(proc.returncode) + \
              'stdout:\n{}\n'.format(stdout_data) + \
              'stderr:\n{}\n'.format(stderr_data)
        assert False, msg
    else:
        if isinstance(stdout_data, bytes):
            stdout_data = stdout_data.decode('utf-8')
        print(stdout_data)

    return

def test_build_mt3dusgs():
    if pymake is None:
        return
    starget = 'MT3D-USGS'
    exe_name = 'mt3dusgs'
    dirname = 'mt3d-usgs_Distribution'
    url = "https://water.usgs.gov/ogw/mt3d-usgs/mt3d-usgs_1.0.zip"

    build_target(starget, exe_name, url, dirname)
    return

def set_compiler(starget):
    fct = fc
    cct = cc
    # parse command line arguments to see if user specified options
    # relative to building the target
    msg = ''
    for idx, arg in enumerate(sys.argv):
        if arg.lower() == '--ifort':
            if len(msg) > 0:
                msg += '\n'
            msg += '{} - '.format(arg.lower()) + \
                   '{} will be built with ifort.'.format(starget)
            fct = 'ifort'
        elif arg.lower() == '--cl':
            if len(msg) > 0:
                msg += '\n'
            msg += '{} - '.format(arg.lower()) + \
                   '{} will be built with cl.'.format(starget)
            cct = 'cl'
        elif arg.lower() == '--clang':
            if len(msg) > 0:
                msg += '\n'
            msg += '{} - '.format(arg.lower()) + \
                   '{} will be built with clang.'.format(starget)
            cct = 'clang'
    if len(msg) > 0:
        print(msg)

    return fct, cct


def build_target(starget, exe_name, url, dirname, srcname='src',
                 replace_function=None, verify=True, keep=True,
                 dble=dbleprec, include_subdirs=False):
    print('Determining if {} needs to be built'.format(starget))
    if platform.system().lower() == 'windows':
        exe_name += '.exe'

    exe_exists = flopy.which(exe_name)
    if exe_exists is not None and keep:
        print('No need to build {}'.format(starget) +
              ' since it exists in the current path')
        return

    fct, cct = set_compiler(starget)

    # set up target
    target = os.path.abspath(os.path.join(bindir, exe_name))

    # get current directory
    cpth = os.getcwd()

    # create temporary path
    dstpth = os.path.join('tempbin')
    print('create...{}'.format(dstpth))
    if not os.path.exists(dstpth):
        os.makedirs(dstpth)
    os.chdir(dstpth)

    # Download the distribution
    pymake.download_and_unzip(url, verify=verify)

    # Set srcdir name
    srcdir = os.path.join(dirname, srcname)

    if replace_function is not None:
        replace_function(srcdir)

    # compile code
    print('compiling...{}'.format(os.path.relpath(target)))
    pymake.main(srcdir, target, fct, cct, makeclean=True,
                expedite=False, dryrun=False, double=dble, debug=False,
                include_subdirs=include_subdirs)

    # change back to original path
    os.chdir(cpth)

    # Clean up downloaded directory
    print('delete...{}'.format(dstpth))
    if os.path.isdir(dstpth):
        shutil.rmtree(dstpth)

    msg = '{} does not exist.'.format(os.path.relpath(target))
    assert os.path.isfile(target), msg

    return

if __name__ == '__main__':
    test_build_mfnwt()
    test_build_mt3dusgs()
