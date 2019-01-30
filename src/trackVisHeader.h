#ifndef __trackVisHeader_h
#define __trackVisHeader_h

struct TRACKVIS_HEADER_V2
{
  char          id_string[6]; // first 5 chars must be "TRACK"
  short int     dim[3];
  float         voxel_size[3];
  float         origin[3];
  short int     n_scalars;
  char          scalar_names[10][20];
  short int     n_properties;
  char          property_names[10][20];
  float         vox_to_ras[4][4];
  char          reserved[444];
  char          voxel_order[4];
  char          pad2[4];
  float         image_orientation_patient[6];
  char          pad1[2];
  unsigned char invert_x;
  unsigned char invert_y;
  unsigned char invert_z;
  unsigned char swap_xy;
  unsigned char swap_yz;
  unsigned char swap_zx;
  int           n_count;
  int           version;
  int           hdr_size;
};

#endif
