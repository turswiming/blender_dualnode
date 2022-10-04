void main()
{
  ivec2 coord_in = ivec2(gl_GlobalInvocationID.xy);
  ivec2 coord_out = coord_in;
  vec4 paint_color = imageLoad(in_paint_img, coord_in);
  paint_color.a = 1.0;
  imageStore(out_img, coord_out, paint_color);
}