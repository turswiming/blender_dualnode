#ifdef USE_WORLD_CLIP_PLANES
#  if defined(GPU_VERTEX_SHADER) || defined(GPU_GEOMETRY_SHADER)

/* When all shaders are builtin shaders are migrated this could be applied directly. */
#    ifdef USE_GPU_SHADER_CREATE_INFO
#      define WorldClipPlanes clipPlanes.world
#    else
uniform vec4 WorldClipPlanes[6];
#    endif

void world_clip_planes_calc_clip_distance(vec3 wpos)
{
  vec4 clip_planes[6] = WorldClipPlanes;
  vec4 pos = vec4(wpos, 1.0);

  gl_ClipDistance[0] = dot(clip_planes[0], pos);
  gl_ClipDistance[1] = dot(clip_planes[1], pos);
  gl_ClipDistance[2] = dot(clip_planes[2], pos);
  gl_ClipDistance[3] = dot(clip_planes[3], pos);
  gl_ClipDistance[4] = dot(clip_planes[4], pos);
  gl_ClipDistance[5] = dot(clip_planes[5], pos);
}

#  endif

#endif
