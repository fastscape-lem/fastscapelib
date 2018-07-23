from _fastscapelib_py import *

versions = get_versions()
__version__ = versions['version']
__git_revision__ = versions['git_hash_full']
del get_versions, versions
