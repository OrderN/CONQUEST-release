/* -----------------------------------------------------------------------------
 * $Id: $
 * -----------------------------------------------------------------------------
 * File cq_ut_conqtour_img.c
 * -----------------------------------------------------------------------------
 *
 * ***** Conquest/utilities/cq_ut_conqtour_img.c *
 *
 * NAME
 *  cq_ut_conqtour_img.c
 * PURPOSE
 *
 * USES
 *
 * AUTHOR
 *  torralba
 * CREATION DATE
 *  Oct 25, 2010
 * MODIFICATION HISTORY
 *
 * *****/

#include "cq_ut_conqtour.h"

void writePpmImage(char *file)
{
  GLubyte *pixels;
  FILE *fimg;
  int i;
  int size;

  fimg = fopen(file, "wb");
  if (fimg == NULL)
    {
      fprintf(stderr, "Problem creating file '%s'\n", file);
      exit(-1);
    }
  size = ctrl.w * ctrl.h * 3;
  pixels = (GLubyte *) malloc((size_t) size * sizeof(GLubyte));
  if (pixels == NULL)
    {
      fprintf(stderr, "Not enough memory for image file buffer\n");
      exit(-1);
    }
  /* Alignment required for 24 bit data (the default is 4, for 32 bits) */
  glPixelStorei(GL_PACK_ALIGNMENT, 1);
  glReadPixels(0, 0, ctrl.w, ctrl.h, GL_RGB, GL_UNSIGNED_BYTE, pixels);
  fprintf(fimg, "P6 %d %d 255\n", ctrl.w, ctrl.h);
  /* The data comes from bottom to top, so flip it */
  for (i = 0; i < ctrl.h; ++i)
    {
      fwrite(&pixels[(ctrl.h - i - 1) * ctrl.w * 3], 1, ctrl.w * 3, fimg);
    }
  fclose(fimg);
  free(pixels);
}
