import gpu
import sys

print('GPU_VENDOR:'+gpu.platform.vendor_get())
print('GPU_RENDERER:'+gpu.platform.renderer_get())
print('GPU_VERSION:'+gpu.platform.version_get())
sys.exit(0)
