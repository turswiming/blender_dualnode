"""
    Prints GPU backend information to the console and exits.
    
    Use this script as `blender --python gpu_info.py`.
    Doesn't work with `--background` parameter as then the GPU backend won't
    be initialized.
"""
import gpu
import sys

print('GPU_VENDOR:' + gpu.platform.vendor_get())
print('GPU_RENDERER:' + gpu.platform.renderer_get())
print('GPU_VERSION:' + gpu.platform.version_get())
print('GPU_DEVICE_TYPE:' + gpu.platform.device_type_get())

sys.exit(0)