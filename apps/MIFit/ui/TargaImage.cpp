#include "TargaImage.h"

#include <stdio.h>

struct TGAHEADER
{
    GLbyte identsize;            // Size of ID field that follows header (0)
    GLbyte colorMapType;         // 0 = None, 1 = paletted
    GLbyte imageType;            // 0 = none, 1 = indexed, 2 = rgb, 3 = grey, +8=rle
    unsigned short colorMapStart;        // First colour map entry
    unsigned short colorMapLength;       // Number of colors
    unsigned char colorMapBits; // bits per palette entry
    unsigned short xstart;               // image x origin
    unsigned short ystart;               // image y origin
    unsigned short width;                // width in pixels
    unsigned short height;               // height in pixels
    GLbyte bits;                 // bits per pixel (8 16, 24, 32)
    GLbyte descriptor;           // image descriptor
};



////////////////////////////////////////////////////////////////////
// Allocate memory and load targa bits. Returns pointer to new buffer,
// height, and width of texture, and the OpenGL format of data.
// Call free() on buffer when finished!
// This only works on pretty vanilla targas... 8, 24, or 32 bit color
// only, no palettes, no RLE encoding.
TargaImage::TargaImage(const std::string &fileName)
{
    FILE *file;    // File pointer
    TGAHEADER tgaHeader;  // TGA file header
    unsigned long imageSize; // Size in bytes of image
    short depth;   // Pixel depth;
    data = NULL;        // Pointer to bits

    // Default/Failed values
    width = 0;
    height = 0;
    format = GL_BGR_EXT;
    components = GL_RGB8;

    // Attempt to open the fil
    file = fopen(fileName.c_str(), "rb");
    if (file == NULL)
    {
        return;
    }

    // Read in header (binary)

    fread(&tgaHeader.identsize, sizeof(GLbyte), 1, file);
    fread(&tgaHeader.colorMapType, sizeof(GLbyte), 1, file);
    fread(&tgaHeader.imageType, sizeof(GLbyte), 1, file);
    fread(&tgaHeader.colorMapStart, sizeof(unsigned short), 1, file);
    fread(&tgaHeader.colorMapLength, sizeof(unsigned short), 1, file);
    fread(&tgaHeader.colorMapBits, sizeof(unsigned char), 1, file);
    fread(&tgaHeader.xstart, sizeof(unsigned short), 1, file);
    fread(&tgaHeader.ystart, sizeof(unsigned short), 1, file);
    fread(&tgaHeader.width, sizeof(unsigned short), 1, file);
    fread(&tgaHeader.height, sizeof(unsigned short), 1, file);
    fread(&tgaHeader.bits, sizeof(GLbyte), 1, file);
    fread(&tgaHeader.descriptor, sizeof(GLbyte), 1, file);

    // Get width, height, and depth of texture
    width = tgaHeader.width;
    height = tgaHeader.height;
    depth = tgaHeader.bits / 8;

    // Put some validity checks here. Very simply, I only understand
    // or care about 8, 24, or 32 bit targa's.
    if (tgaHeader.bits != 8 && tgaHeader.bits != 24 && tgaHeader.bits != 32)
    {
        return;
    }

    // Calculate size of image buffer
    imageSize = tgaHeader.width * tgaHeader.height * depth;

    // Allocate memory and check for success
    data = new GLubyte[imageSize];
    if (data == NULL)
    {
        return;
    }

    // Read in the bits
    // Check for read error. This should catch RLE or other
    // weird formats that I don't want to recognize
    if (fread(data, imageSize, 1, file) != 1)
    {
        delete[] data;
        data = NULL;
        return;
    }

    // Set OpenGL format expected
    switch (depth)
    {
    case 3:   // Most likely case
        format = GL_BGR_EXT;
        components = GL_RGB8;
        break;
    case 4:
        format = GL_BGRA_EXT;
        components = GL_RGBA8;
        break;
    case 1:
        format = GL_LUMINANCE;
        components = GL_LUMINANCE8;
        break;
    }


    // Done with File
    fclose(file);
}

TargaImage::~TargaImage()
{
    if (data != NULL)
    {
        delete[] data;
    }
}

