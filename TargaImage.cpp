///////////////////////////////////////////////////////////////////////////////
//
//      TargaImage.cpp                          Author:     Stephen Chenney
//                                              Modified:   Eric McDaniel
//                                              Date:       Fall 2004
//                                              Modified:   Feng Liu
//                                              Date:       Winter 2011
//                                              Why:        Change the library file 
//      Implementation of TargaImage methods.  You must implement the image
//  modification functions.
//
///////////////////////////////////////////////////////////////////////////////

#include "Globals.h"
#include "TargaImage.h"
#include "libtarga.h"
#include <stdlib.h>
#include <assert.h>
#include <memory.h>
#include <math.h>
#include <iostream>
#include <sstream>
#include <vector>
#include <algorithm>

using namespace std;

// constants
const int           RED             = 0;                // red channel
const int           GREEN           = 1;                // green channel
const int           BLUE            = 2;                // blue channel
const unsigned char BACKGROUND[3]   = { 0, 0, 0 };      // background color


// Computes n choose s, efficiently
double Binomial(int n, int s)
{
    double        res;

    res = 1;
    for (int i = 1 ; i <= s ; i++)
        res = (n - i + 1) * res / i ;

    return res;
}// Binomial


///////////////////////////////////////////////////////////////////////////////
//
//      Constructor.  Initialize member variables.
//
///////////////////////////////////////////////////////////////////////////////
TargaImage::TargaImage() : width(0), height(0), data(NULL)
{}// TargaImage

///////////////////////////////////////////////////////////////////////////////
//
//      Constructor.  Initialize member variables.
//
///////////////////////////////////////////////////////////////////////////////
TargaImage::TargaImage(int w, int h) : width(w), height(h)
{
   data = new unsigned char[width * height * 4];
   ClearToBlack();
}// TargaImage



///////////////////////////////////////////////////////////////////////////////
//
//      Constructor.  Initialize member variables to values given.
//
///////////////////////////////////////////////////////////////////////////////
TargaImage::TargaImage(int w, int h, unsigned char *d)
{
    int i;

    width = w;
    height = h;
    data = new unsigned char[width * height * 4];

    for (i = 0; i < width * height * 4; i++)
	    data[i] = d[i];
}// TargaImage

///////////////////////////////////////////////////////////////////////////////
//
//      Copy Constructor.  Initialize member to that of input
//
///////////////////////////////////////////////////////////////////////////////
TargaImage::TargaImage(const TargaImage& image) 
{
   width = image.width;
   height = image.height;
   data = NULL; 
   if (image.data != NULL) {
      data = new unsigned char[width * height * 4];
      memcpy(data, image.data, sizeof(unsigned char) * width * height * 4);
   }
}


///////////////////////////////////////////////////////////////////////////////
//
//      Destructor.  Free image memory.
//
///////////////////////////////////////////////////////////////////////////////
TargaImage::~TargaImage()
{
    if (data)
        delete[] data;
}// ~TargaImage


///////////////////////////////////////////////////////////////////////////////
//
//      Converts an image to RGB form, and returns the rgb pixel data - 24 
//  bits per pixel. The returned space should be deleted when no longer 
//  required.
//
///////////////////////////////////////////////////////////////////////////////
unsigned char* TargaImage::To_RGB(void)
{
    unsigned char   *rgb = new unsigned char[width * height * 3];
    int		    i, j;

    if (! data)
	    return NULL;

    // Divide out the alpha
    for (i = 0 ; i < height ; i++)
    {
	    int in_offset = i * width * 4;
	    int out_offset = i * width * 3;

	    for (j = 0 ; j < width ; j++)
        {
	        RGBA_To_RGB(data + (in_offset + j*4), rgb + (out_offset + j*3));
	    }
    }

    return rgb;
}// TargaImage


///////////////////////////////////////////////////////////////////////////////
//
//      Save the image to a targa file. Returns 1 on success, 0 on failure.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Save_Image(const char *filename)
{
    TargaImage	*out_image = Reverse_Rows();

    if (! out_image)
	    return false;

    if (!tga_write_raw(filename, width, height, out_image->data, TGA_TRUECOLOR_32))
    {
	    cout << "TGA Save Error: %s\n", tga_error_string(tga_get_last_error());
	    return false;
    }

    delete out_image;

    return true;
}// Save_Image


///////////////////////////////////////////////////////////////////////////////
//
//      Load a targa image from a file.  Return a new TargaImage object which 
//  must be deleted by caller.  Return NULL on failure.
//
///////////////////////////////////////////////////////////////////////////////
TargaImage* TargaImage::Load_Image(char *filename)
{
    unsigned char   *temp_data;
    TargaImage	    *temp_image;
    TargaImage	    *result;
    int		        width, height;

    if (!filename)
    {
        cout << "No filename given." << endl;
        return NULL;
    }// if

    temp_data = (unsigned char*)tga_load(filename, &width, &height, TGA_TRUECOLOR_32);
    if (!temp_data)
    {
        cout << "TGA Error: %s\n", tga_error_string(tga_get_last_error());
	    width = height = 0;
	    return NULL;
    }
    temp_image = new TargaImage(width, height, temp_data);
    free(temp_data);

    result = temp_image->Reverse_Rows();

    delete temp_image;

    return result;
}// Load_Image


///////////////////////////////////////////////////////////////////////////////
//
//      Convert image to grayscale.  Red, green, and blue channels should all 
//  contain grayscale value.  Alpha channel shoould be left unchanged.  Return
//  success of operation.
/*
    Use the formula I = 0.299r + 0.587g + 0.114b to convert color images to grayscale. 
    This will be a key pre-requisite for many other operations. This operation should 
    not affect alpha in any way.
*/
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::To_Grayscale()
{
    if (data == NULL)
        return false;

    int total_pixels = width * height;
    int size = total_pixels * 4; // RGBA

    unsigned char red = 0;
    unsigned char blue = 0;
    unsigned char green = 0;
    unsigned char intensity = 0;

    for (int i = 0; i < (size - 4); i = i + 4)
    {
        red = data[i];
        green = data[i + 1];
        blue = data[i + 2];

        // Method from project specification: Weighted average
        /* I = 0.299r + 0.587g + 0.114b */
        intensity = (red * 0.299) + (green * 0.587) + (blue * 0.114);

        // R = G = B = gray
        data[i] = intensity;
        data[i + 1] = intensity;
        data[i + 2] = intensity;
    }
    return true;
}// To_Grayscale


///////////////////////////////////////////////////////////////////////////////
//
//  Convert the image to an 8 bit image using uniform quantization.  Return 
//  success of operation.
//
/*
	All of these operations assume that the current image has 24 bits of color 
	information. They should still produce 24 bit images, but there should only 
	be 256 different colors in the resulting image (so the image could be stored 
	as an 8 bit indexed color image). Don't be concerned with what happens if you 
	run these operations on something that is already quantized. These operations 
	should not affect alpha - we will only test them on images with 
	alpha = 1 (fully opaque images).

    Use the uniform quantization algorithm to convert the current image from a 24 
    bit color image to an 8 bit color image. Use 4 levels of blue, 8 levels of red, 
    and 8 levels of green in the quantized image.

*/
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Quant_Uniform()
{
    if(data == NULL)
		return false;

    int total_pixels = width * height;
    int size = total_pixels * 4;

    // These are used to tally all color totals 
    unsigned char r_value_count[256]{}; 
    unsigned char g_value_count[256]{}; 
    unsigned char b_value_count[256]{}; 

    // New color spaces (256 colors)
    unsigned char r_space[8]{};
    unsigned char g_space[8]{};
    unsigned char b_space[4]{};

    const int rg_offset = 32;
    const int b_offset = 64;

    // counts number of colors for each channel
    for (int i = 0; i < (size - 4); i = i + 4)
    {
        ++r_value_count[data[i]];
        ++g_value_count[data[i + 1]];
        ++b_value_count[data[i + 2]];
    }

    // gets the average colors for rg space
    for (int i = 0, j = 0; i < 256; i = i + rg_offset, j++)
    {
        int r_sum = 0;
        int r_count = 0;

        int g_sum = 0;
        int g_count = 0;

        unsigned char r_average = 0;
        unsigned char g_average = 0;

        // counting and summing each range for rg_space
        for (int value = i; value < (i + rg_offset); value++)
        {
            r_sum += (value * r_value_count[value]);
            r_count += r_value_count[value];

            g_sum += (value * g_value_count[value]);
            g_count += g_value_count[value];
        }

        r_average = (unsigned char) (r_sum / r_count);
        g_average = (unsigned char) (g_sum / g_count);

        // save local average to build up color space
        r_space[j] = r_average;
        g_space[j] = g_average;
    }

    // gets the average colors for b space (new offest)
    for (int i = 0, j = 0; i < 256; i = i + b_offset, j++)
    {
        int b_sum = 0;
        int b_count = 0;

        unsigned char b_average = 0;

        // counting and summing each range for b_space
        for (int value = i; value < (i + b_offset); value++)
        {
            b_sum += (value * b_value_count[value]);
            b_count += b_value_count[value];
        }

        b_average = (unsigned char) (b_sum / b_count);

        // save local average to build up color space
        b_space[j] = b_average;
    }

    // Applies the uniform quantization 
    int r = 0;
    int g = 0;
    int b = 0;
    for (int i = 0; i < (size - 4); i = i + 4)
    {
        r = (int) (data[i] / rg_offset);
        g = (int) (data[i + 1] / rg_offset);
        b = (int) (data[i + 2] / b_offset);

        data[i] = r_space[r];
        data[i + 1] = g_space[g];
        data[i + 2] = b_space[b];
    }

    // Uncomment to print new color spaces to conosole.
		std::cout << "New R space: ";
		LogColorSpace(r_space, 8);

		std::cout << "New G space: ";
		LogColorSpace(g_space, 8);

		std::cout << "New B space: ";
		LogColorSpace(b_space, 4);

    return true;
}// Quant_Uniform

// Prints space array typecasted as ints to the console. 
void TargaImage::LogColorSpace(const unsigned char* space, const int size) const
{
    std::cout << "[";
    for (int i = 0; i < size; i++)
    {
        std::cout << (int) space[i];
        if (i != (size - 1))
            std::cout << ", ";
        else
            std::cout << "]" << std::endl;
    }
}


///////////////////////////////////////////////////////////////////////////////
//
//      Convert the image to an 8 bit image using populosity quantization.  
//  Return success of operation.
/*
        Use the populosity algorithm to convert the current 24 bit color image 
        to an 8 bit color image. Before building the color usage histogram, do 
        a uniform quantization step down to 32 levels of each primary. This gives 
        32 x 32 x 32 = 32768 possible colors. Then find the 256 most popular colors, 
        then map the original colors onto their closest chosen color. To find the 
        closest color, use the euclidean (L2) distance in RGB space. If (r1,g1,b1) 
        and (r2,g2,b2) are the colors, use sqrt((r1-r2)^2 + (g1-g2)^2 + (b1-b2)^2) 
        suitably converted into C++ code.

    // 1. Do a uniform-quant step down to 32 levels of each primary color. 
    // 2. Find 256 most popular colors
    // 3. Map original colors onto closest chosen color
*/
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Quant_Populosity()
{
    if(data == NULL)
		return false;

    int total_pixels = width * height;
    int size = total_pixels * 4;

    // Each index represents a color value, and each element represents that
    // color's count in the current image.
    unsigned char r_value_count[256]{}; 
    unsigned char g_value_count[256]{}; 
    unsigned char b_value_count[256]{}; 

    // 32 levels for each primary color
    unsigned char r_space[32]{};
    unsigned char g_space[32]{};
    unsigned char b_space[32]{};

    // sliding window for taking the average over value_count arrays
    const int rgb_offset = 8;

    // counts number of colors for each channel
    for (int i = 0; i < (size - 4); i = i + 4)
    {
        ++r_value_count[data[i]];
        ++g_value_count[data[i + 1]];
        ++b_value_count[data[i + 2]];
        // i + 3 == alpha channel, so we don't touch 
    }

    // Traverses the 256 rgb_value_count arrays and takes a sliding
    // average across with a radius of rbg_offset.
    for (int i = 0, j = 0; i < 256; i = i + rgb_offset, j++)
    {
        int r_sum = 0;
        int r_count = 0;

        int g_sum = 0;
        int g_count = 0;

        int b_sum = 0;
        int b_count = 0;

        unsigned char r_average = 0;
        unsigned char g_average = 0;
        unsigned char b_average = 0;

        // counts and sums each channel within sliding window
        for (int value = i; value < (i + rgb_offset); value++)
        {
            r_sum += (value * r_value_count[value]);
            r_count += r_value_count[value];

            g_sum += (value * g_value_count[value]);
            g_count += g_value_count[value];

            b_sum += (value * b_value_count[value]);
            b_count += b_value_count[value];
        }

        // calculate the average from traversing window
        r_average = (unsigned char) (r_sum / r_count);
        g_average = (unsigned char) (g_sum / g_count);
        b_average = (unsigned char) (b_sum / b_count);

        // save local average to build up color space
        r_space[j] = r_average;
        g_space[j] = g_average;
        b_space[j] = b_average;
    }

    // Applies the uniform quantization 
    int r = 0;
    int g = 0;
    int b = 0;
    for (int i = 0; i < (size - 4); i = i + 4)
    {
        // sine the rgb_space is quantizied to 32 spaces if we divide
        // the original image data by the rgb_offset (32) then we get
        // the index for the appropriate color to reassign.
        r = (int) (data[i] / rgb_offset);
        g = (int) (data[i + 1] / rgb_offset);
        b = (int) (data[i + 2] / rgb_offset);

        data[i] = r_space[r];
        data[i + 1] = g_space[g];
        data[i + 2] = b_space[b];
        // data[i+3] is the alpha channel
    }

    // Uncomment to print new color spaces to conosole.
    /*
    std::cout << "New R space: ";
	LogColorSpace(r_space, 32);

	std::cout << "New G space: ";
	LogColorSpace(g_space, 32);

	std::cout << "New B space: ";
	LogColorSpace(b_space, 32);
    */


    // Now we have quantizied to 32 levels of RGB!

    // Now we need to find the 256 most popular colors...


    return true;
}// Quant_Populosity


///////////////////////////////////////////////////////////////////////////////
//
//      Dither the image using a threshold of 1/2.  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Dither_Threshold()
{
    ClearToBlack();
    return false;
}// Dither_Threshold


///////////////////////////////////////////////////////////////////////////////
//
//      Dither image using random dithering.  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Dither_Random()
{
    ClearToBlack();
    return false;
}// Dither_Random


///////////////////////////////////////////////////////////////////////////////
//
//      Perform Floyd-Steinberg dithering on the image.  Return success of 
//  operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Dither_FS()
{
    ClearToBlack();
    return false;
}// Dither_FS


///////////////////////////////////////////////////////////////////////////////
//
//      Dither the image while conserving the average brightness.  Return 
//  success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Dither_Bright()
{
    ClearToBlack();
    return false;
}// Dither_Bright


///////////////////////////////////////////////////////////////////////////////
//
//      Perform clustered differing of the image.  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Dither_Cluster()
{
    ClearToBlack();
    return false;
}// Dither_Cluster


///////////////////////////////////////////////////////////////////////////////
//
//  Convert the image to an 8 bit image using Floyd-Steinberg dithering over
//  a uniform quantization - the same quantization as in Quant_Uniform.
//  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Dither_Color()
{
    ClearToBlack();
    return false;
}// Dither_Color


///////////////////////////////////////////////////////////////////////////////
//
//      Composite the current image over the given image.  Return success of 
//  operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Comp_Over(TargaImage* pImage)
{
    if (width != pImage->width || height != pImage->height)
    {
        cout <<  "Comp_Over: Images not the same size\n";
        return false;
    }

    ClearToBlack();
    return false;
}// Comp_Over


///////////////////////////////////////////////////////////////////////////////
//
//      Composite this image "in" the given image.  See lecture notes for 
//  details.  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Comp_In(TargaImage* pImage)
{
    if (width != pImage->width || height != pImage->height)
    {
        cout << "Comp_In: Images not the same size\n";
        return false;
    }

    ClearToBlack();
    return false;
}// Comp_In


///////////////////////////////////////////////////////////////////////////////
//
//      Composite this image "out" the given image.  See lecture notes for 
//  details.  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Comp_Out(TargaImage* pImage)
{
    if (width != pImage->width || height != pImage->height)
    {
        cout << "Comp_Out: Images not the same size\n";
        return false;
    }

    ClearToBlack();
    return false;
}// Comp_Out


///////////////////////////////////////////////////////////////////////////////
//
//      Composite current image "atop" given image.  Return success of 
//  operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Comp_Atop(TargaImage* pImage)
{
    if (width != pImage->width || height != pImage->height)
    {
        cout << "Comp_Atop: Images not the same size\n";
        return false;
    }

    ClearToBlack();
    return false;
}// Comp_Atop


///////////////////////////////////////////////////////////////////////////////
//
//      Composite this image with given image using exclusive or (XOR).  Return
//  success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Comp_Xor(TargaImage* pImage)
{
    if (width != pImage->width || height != pImage->height)
    {
        cout << "Comp_Xor: Images not the same size\n";
        return false;
    }

    ClearToBlack();
    return false;
}// Comp_Xor


///////////////////////////////////////////////////////////////////////////////
//
//      Calculate the difference bewteen this imag and the given one.  Image 
//  dimensions must be equal.  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Difference(TargaImage* pImage)
{
    if (!pImage)
        return false;

    if (width != pImage->width || height != pImage->height)
    {
        cout << "Difference: Images not the same size\n";
        return false;
    }// if

    for (int i = 0 ; i < width * height * 4 ; i += 4)
    {
        unsigned char        rgb1[3];
        unsigned char        rgb2[3];

        RGBA_To_RGB(data + i, rgb1);
        RGBA_To_RGB(pImage->data + i, rgb2);

        data[i] = abs(rgb1[0] - rgb2[0]);
        data[i+1] = abs(rgb1[1] - rgb2[1]);
        data[i+2] = abs(rgb1[2] - rgb2[2]);
        data[i+3] = 255;
    }

    return true;
}// Difference


///////////////////////////////////////////////////////////////////////////////
//
//      Perform 5x5 box filter on this image.  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Filter_Box()
{
    ClearToBlack();
    return false;
}// Filter_Box


///////////////////////////////////////////////////////////////////////////////
//
//      Perform 5x5 Bartlett filter on this image.  Return success of 
//  operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Filter_Bartlett()
{
    ClearToBlack();
    return false;
}// Filter_Bartlett


///////////////////////////////////////////////////////////////////////////////
//
//      Perform 5x5 Gaussian filter on this image.  Return success of 
//  operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Filter_Gaussian()
{
    ClearToBlack();
    return false;
}// Filter_Gaussian

///////////////////////////////////////////////////////////////////////////////
//
//      Perform NxN Gaussian filter on this image.  Return success of 
//  operation.
//
///////////////////////////////////////////////////////////////////////////////

bool TargaImage::Filter_Gaussian_N( unsigned int N )
{
    ClearToBlack();
   return false;
}// Filter_Gaussian_N


///////////////////////////////////////////////////////////////////////////////
//
//      Perform 5x5 edge detect (high pass) filter on this image.  Return 
//  success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Filter_Edge()
{
    ClearToBlack();
    return false;
}// Filter_Edge


///////////////////////////////////////////////////////////////////////////////
//
//      Perform a 5x5 enhancement filter to this image.  Return success of 
//  operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Filter_Enhance()
{
    ClearToBlack();
    return false;
}// Filter_Enhance


///////////////////////////////////////////////////////////////////////////////
//
//      Run simplified version of Hertzmann's painterly image filter.
//      You probably will want to use the Draw_Stroke funciton and the
//      Stroke class to help.
// Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::NPR_Paint()
{
    ClearToBlack();
    return false;
}



///////////////////////////////////////////////////////////////////////////////
//
//      Halve the dimensions of this image.  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Half_Size()
{
    ClearToBlack();
    return false;
}// Half_Size


///////////////////////////////////////////////////////////////////////////////
//
//      Double the dimensions of this image.  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Double_Size()
{
    ClearToBlack();
    return false;
}// Double_Size


///////////////////////////////////////////////////////////////////////////////
//
//      Scale the image dimensions by the given factor.  The given factor is 
//  assumed to be greater than one.  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Resize(float scale)
{
    ClearToBlack();
    return false;
}// Resize


//////////////////////////////////////////////////////////////////////////////
//
//      Rotate the image clockwise by the given angle.  Do not resize the 
//  image.  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Rotate(float angleDegrees)
{
    ClearToBlack();
    return false;
}// Rotate


//////////////////////////////////////////////////////////////////////////////
//
//      Given a single RGBA pixel return, via the second argument, the RGB
//      equivalent composited with a black background.
//
///////////////////////////////////////////////////////////////////////////////
void TargaImage::RGBA_To_RGB(unsigned char *rgba, unsigned char *rgb)
{
    const unsigned char	BACKGROUND[3] = { 0, 0, 0 };

    unsigned char  alpha = rgba[3];

    if (alpha == 0)
    {
        rgb[0] = BACKGROUND[0];
        rgb[1] = BACKGROUND[1];
        rgb[2] = BACKGROUND[2];
    }
    else
    {
	    float	alpha_scale = (float)255 / (float)alpha;
	    int	val;
	    int	i;

	    for (i = 0 ; i < 3 ; i++)
	    {
	        val = (int)floor(rgba[i] * alpha_scale);
	        if (val < 0)
		    rgb[i] = 0;
	        else if (val > 255)
		    rgb[i] = 255;
	        else
		    rgb[i] = val;
	    }
    }
}// RGA_To_RGB


///////////////////////////////////////////////////////////////////////////////
//
//      Copy this into a new image, reversing the rows as it goes. A pointer
//  to the new image is returned.
//
///////////////////////////////////////////////////////////////////////////////
TargaImage* TargaImage::Reverse_Rows(void)
{
    unsigned char   *dest = new unsigned char[width * height * 4];
    TargaImage	    *result;
    int 	        i, j;

    if (! data)
    	return NULL;

    for (i = 0 ; i < height ; i++)
    {
	    int in_offset = (height - i - 1) * width * 4;
	    int out_offset = i * width * 4;

	    for (j = 0 ; j < width ; j++)
        {
	        dest[out_offset + j * 4] = data[in_offset + j * 4];
	        dest[out_offset + j * 4 + 1] = data[in_offset + j * 4 + 1];
	        dest[out_offset + j * 4 + 2] = data[in_offset + j * 4 + 2];
	        dest[out_offset + j * 4 + 3] = data[in_offset + j * 4 + 3];
        }
    }

    result = new TargaImage(width, height, dest);
    delete[] dest;
    return result;
}// Reverse_Rows


///////////////////////////////////////////////////////////////////////////////
//
//      Clear the image to all black.
//
///////////////////////////////////////////////////////////////////////////////
void TargaImage::ClearToBlack()
{
    memset(data, 0, width * height * 4);
}// ClearToBlack


///////////////////////////////////////////////////////////////////////////////
//
//      Helper function for the painterly filter; paint a stroke at
// the given location
//
///////////////////////////////////////////////////////////////////////////////
void TargaImage::Paint_Stroke(const Stroke& s) {
   int radius_squared = (int)s.radius * (int)s.radius;
   for (int x_off = -((int)s.radius); x_off <= (int)s.radius; x_off++) {
      for (int y_off = -((int)s.radius); y_off <= (int)s.radius; y_off++) {
         int x_loc = (int)s.x + x_off;
         int y_loc = (int)s.y + y_off;
         // are we inside the circle, and inside the image?
         if ((x_loc >= 0 && x_loc < width && y_loc >= 0 && y_loc < height)) {
            int dist_squared = x_off * x_off + y_off * y_off;
            if (dist_squared <= radius_squared) {
               data[(y_loc * width + x_loc) * 4 + 0] = s.r;
               data[(y_loc * width + x_loc) * 4 + 1] = s.g;
               data[(y_loc * width + x_loc) * 4 + 2] = s.b;
               data[(y_loc * width + x_loc) * 4 + 3] = s.a;
            } else if (dist_squared == radius_squared + 1) {
               data[(y_loc * width + x_loc) * 4 + 0] = 
                  (data[(y_loc * width + x_loc) * 4 + 0] + s.r) / 2;
               data[(y_loc * width + x_loc) * 4 + 1] = 
                  (data[(y_loc * width + x_loc) * 4 + 1] + s.g) / 2;
               data[(y_loc * width + x_loc) * 4 + 2] = 
                  (data[(y_loc * width + x_loc) * 4 + 2] + s.b) / 2;
               data[(y_loc * width + x_loc) * 4 + 3] = 
                  (data[(y_loc * width + x_loc) * 4 + 3] + s.a) / 2;
            }
         }
      }
   }
}



///////////////////////////////////////////////////////////////////////////////
//
//      Build a Stroke
//
///////////////////////////////////////////////////////////////////////////////
Stroke::Stroke() {}

///////////////////////////////////////////////////////////////////////////////
//
//      Build a Stroke
//
///////////////////////////////////////////////////////////////////////////////
Stroke::Stroke(unsigned int iradius, unsigned int ix, unsigned int iy,
               unsigned char ir, unsigned char ig, unsigned char ib, unsigned char ia) :
   radius(iradius),x(ix),y(iy),r(ir),g(ig),b(ib),a(ia)
{
}

