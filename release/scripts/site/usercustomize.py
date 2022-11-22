# SPDX-License-Identifier: GPL-2.0-or-later
# Make shared libraries needed by modules available in standalone Python binary.

import sys
import os

if sys.platform == 'win32':
    exe_dir, exe_file = os.path.split(sys.executable)
    if exe_file.startswith('python'):
        blender_dir = os.path.abspath(os.path.join(exe_dir, '..', '..', '..','blender.shared'))
        os.add_dll_directory(blender_dir)
        # OIIO will by default add all paths from the path variable to add_dll_directory
        # problem there is that those folders will be searched before ours and versions of
        # some dlls may be found that are not blenders and may not even be the right version
        # causing compatibility issues.
        os.environ["OIIO_LOAD_DLLS_FROM_PATH"] = "0"

        import_paths = os.getenv('PXR_USD_WINDOWS_DLL_PATH')
        if import_paths is None:
            os.environ["PXR_USD_WINDOWS_DLL_PATH"] = blender_dir

        materialx_libs_dir = os.path.abspath(os.path.join(exe_dir, '..', '..', 'datafiles', 'materialx', 'libraries'))
        materialx_libs_env = os.getenv('MATERIALX_SEARCH_PATH')
        if materialx_libs_env is None:
            os.environ["MATERIALX_SEARCH_PATH"] = materialx_libs_dir
