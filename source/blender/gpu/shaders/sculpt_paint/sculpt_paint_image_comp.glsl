void main()
{
  PackedPixelRow row = pixel_row_buf[gl_GlobalInvocationID.x + pixel_row_offset];
  ivec2 image_coord = PIXEL_ROW_START_IMAGE_COORD(row);

  uint row_len = PIXEL_ROW_LEN(row);

  for (int i = 0; i < row_len; i++) {

    imageStore(out_img, image_coord + ivec2(i, 0), float4(1.0, 0.0, 1.0, 1.0));
  }
}