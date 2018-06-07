#ifndef IMAGE_WRITER_H_1
#define IMAGE_WRITER_H_1
#ifdef __APPLE__
#include <GLUT/glut.h>
#include <OpenGL/gl.h>
#include <OpenGL/glu.h>
#else
#include <GL/glut.h>
#endif
#ifdef WIN32
#include <windows.h>
#endif
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

//Image class 
//This class can be used to write your final result to an image. 
//You can open the image using a PPM viewer.

//YOU CAN IGNORE THIS CODE!
 class RGBValue
{
	public:
	RGBValue(float rI=0, float gI=0, float bI=0)
	: r(rI)
	, g(gI)
	, b(bI)
	{
		if (r>1)
			r=1.0;
		if (g>1)
			g=1.0;
		if (b>1)
			b=1.0;

		if (r<0)
			r=0.0;
		if (g<0)
			g=0.0;
		if (b<0)
			b=0.0;
	};
	
	float operator[](int i) const
	{
		switch(i)
		{
			case 0:
				return r;
			case 1:
				return g;
			case 2:
				return b;
			default: 
				return r;
		}
	}
	float & operator[](int i)
	{
		switch(i)
		{
			case 0:
				return r;
			case 1:
				return g;
			case 2:
				return b;
			default: 
				return r;
		}
	}
	float r, b,g;
};


class Image
{
	public:
	Image(int width, int height)
	: _width(width)
	, _height(height)
	{
		_image.resize(3*_width*_height);
	}
	void setPixel(int i, int j, const RGBValue & rgb)
	{
		_image[3*(_width*j+i)]=rgb[0];
		_image[3*(_width*j+i)+1]=rgb[1];
		_image[3*(_width*j+i)+2]=rgb[2];
		
	}
	std::vector<float> _image;
	int _width;
	int _height;

	bool writeImage(const char * filename);	
};

bool Image::writeImage(const char * filename)
{
    std::string realFilename(filename);
    // Prepend data file directory, so files are in-source.
    realFilename = std::string(OUTPUT_DIR) + '\\' + realFilename;
    // we replace the \ by /
    for (unsigned int i = 0; i<realFilename.length(); ++i)
    {
        if (realFilename[i] == '\\')
            realFilename[i] = '/';
    }

    // Convert everything to unsigned char for image data
    std::vector<unsigned char> imageC(_image.size());
    for (unsigned int i = 0; i < _image.size(); ++i)
    {
        imageC[i] = (unsigned char)(_image[i] * 255.0f);
    }

    // Write and return success/failure
    return stbi_write_bmp(realFilename.c_str(), _width, _height, 3, &imageC[0]);
}

#endif
