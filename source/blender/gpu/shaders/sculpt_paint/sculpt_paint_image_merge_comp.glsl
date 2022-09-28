void main()
{
  ivec2 coord_in = ivec2(gl_GlobalInvocationID.xy);
  ivec2 coord_out = coord_in;
  vec4 paint_color = imageLoad(in_paint_img, coord_in);
  vec4 canvas_color = imageLoad(out_img, coord_out);
  vec4 out_color = mix(canvas_color, paint_color, 0.001);
  imageStore(out_img, coord_out, out_color);
}