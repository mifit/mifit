#ifndef MIFIT_UI_TARGAIMAGE_H_
#define MIFIT_UI_TARGAIMAGE_H_

#include <opengl/OpenGL.h>
#include <string>

class TargaImage
{

    GLint width;
    GLint height;
    GLint components;
    GLenum format;
    GLubyte *data;

public:

    TargaImage(const std::string &fileName);
    ~TargaImage();

    GLint getWidth();
    GLint getHeight();
    GLint getComponents();
    GLenum getFormat();
    GLubyte *getData();

};

inline GLint TargaImage::getWidth()
{
    return width;
}

inline GLint TargaImage::getHeight()
{
    return height;
}

inline GLint TargaImage::getComponents()
{
    return components;
}

inline GLenum TargaImage::getFormat()
{
    return format;
}

inline GLubyte*TargaImage::getData()
{
    return data;
}

#endif /*MIFIT_UI_TARGAIMAGE_H_*/
