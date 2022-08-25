# SPDX-License-Identifier: GPL-2.0-or-later
# Make shared libraries needed by modules available in standalone Python binary.

import sys
import os

if sys.platform == 'win32':
    exe_dir, exe_file = os.path.split(sys.executable)
    if exe_file.startswith('python'):
        blender_dir = os.path.abspath(os.path.join(exe_dir, '..', '..', '..','blender.shared'))
        os.add_dll_directory(blender_dir)
        import_paths = os.getenv('PXR_USD_WINDOWS_DLL_PATH')
        if import_paths is None:
            os.environ["PXR_USD_WINDOWS_DLL_PATH"] = blender_dir
