/*-----[--.----+----.----+----.-----------------------------------------*/

/*      * * F E A P * * A Finite Element Analysis Program               */

/*....  Copyright (c) 1984-2014: Regents of the University of California*/
/*                               All rights reserved                    */

/*----[--.----+----.----+----.-----------------------------------------]
 *    Modification log                                Date (dd-mm-year)
 *      Original version                                    27-04-2011
 */
/*-----[--.----+----.----+----.-----------------------------------------*/
/*     Purpose: JPEG file writer                                        */

/*     Inputs:                                                          */

/*     Outputs:                                                         */
/*-----[--.----+----.----+----.-----------------------------------------*/


#include <stdio.h>                /* Unix standard I/O definitions */
#include <stdlib.h>               /* Use "malloc", "calloc" and "free" */

#include <X11/Xlib.h>
#include <X11/Xutil.h>

/* #include <sys/select.h>  Use this include instead of next 3
                            on POSIX systems                     */
#include <sys/types.h> 
#include <sys/time.h> 
#include <unistd.h> 

#include <jpeglib.h>

#include "digwin.h"


/* Global variables for creating movie file names */
static int animi_jpg = 1;


/* Output a copy of the current graphics window to a JPEG file */
void jpgd(DIGWin* current_dw)

{
  int i,j,px,k,r_stride,rv,gv,bv;

  /* Output file name format FeapXXXXX.jpg */
  char fname[14];           

  XImage *myimage;

  JSAMPLE *jpgdata;
  JSAMPROW rowp[1];
  struct jpeg_compress_struct cinfo;
  struct jpeg_error_mgr jerr;

  FILE *fp;

  /* Create File Name: N.B. Starts looking from Feap00001.jpg 
   * if first time through
   */
  if(animi_jpg == 1) {
    sprintf(fname,"Feap%05d.jpg", animi_jpg++); 
    while( !access(fname,F_OK) ) sprintf(fname,"Feap%05d.jpg", animi_jpg++); 
  }
  else sprintf(fname,"Feap%05d.jpg",animi_jpg++); 

  /* Open File */
  printf("Opening JPEG File %s\n",fname);
  fp = fopen(fname,"w"); 

  /* Modify name for header write */
  fname[9]='_';

  /* Get image from pixmap */
  myimage = XGetImage(current_dw->xdisplay, current_dw->svimage, 
                      0,0, current_dw->xwa.width, current_dw->xwa.height,  
                      AllPlanes , ZPixmap);

  /* Allocate necessary room for RGB data without Alpha */
  jpgdata = (unsigned char *)calloc(current_dw->xwa.width  * 
                                    current_dw->xwa.height * 3,
                                    sizeof(JSAMPLE));

  /* Grab Pixels, remove leading 8 bits, output RGB bits */
  k = 0;
  for(j=0;j<current_dw->xwa.height;j++) {
    for(i=0;i<current_dw->xwa.width;i++) {
	px = (int)XGetPixel(myimage,i,j);

        /* Swap black and white */
        if (px == 0x00ffffff) px = 0x00000000;
        else if (px == 0x00000000) px = 0x00ffffff;

        /* Extract RBG values */
        rv = (px & 0x00ff0000) >> 16;
        gv = (px & 0x0000ff00) >>  8;
        bv = (px & 0x000000ff)      ;

        /* Write out RBG */
        jpgdata[k]=(JSAMPLE)rv;
        jpgdata[k+1]=(JSAMPLE)gv;
        jpgdata[k+2]=(JSAMPLE)bv;
        k = k + 3;
    }
  } 

  cinfo.err = jpeg_std_error(&jerr);
  jpeg_create_compress(&cinfo);
  jpeg_stdio_dest(&cinfo, fp);

  cinfo.image_width = current_dw->xwa.width;
  cinfo.image_height = current_dw->xwa.height;

  cinfo.input_components = 3;           /* # of color components per pixel */
  cinfo.in_color_space = JCS_RGB;       /* colorspace of input image */

  jpeg_set_defaults(&cinfo);

  /* Set quality factor to 100 for lowest compression losses */
  jpeg_set_quality(&cinfo, 100, TRUE);

  /* Set sampling for each color component to 1 for good definition */
  cinfo.comp_info[0].h_samp_factor = 1;
  cinfo.comp_info[0].v_samp_factor = 1;
  cinfo.comp_info[1].h_samp_factor = 1;
  cinfo.comp_info[1].h_samp_factor = 1;
  cinfo.comp_info[2].v_samp_factor = 1;
  cinfo.comp_info[2].v_samp_factor = 1;

  jpeg_start_compress(&cinfo, TRUE);

  r_stride = current_dw->xwa.width * 3; /* JSAMPLEs per row in image_buffer */

  while (cinfo.next_scanline < cinfo.image_height) {
    rowp[0] = & jpgdata[cinfo.next_scanline * r_stride];
    (void) jpeg_write_scanlines(&cinfo, rowp, 1);
  }

  jpeg_finish_compress(&cinfo);
  jpeg_destroy_compress(&cinfo);

  free(jpgdata);
  fclose(fp);
}

