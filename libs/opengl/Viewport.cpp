#include <opengl/Viewport.h>

#include <opengl/OpenGL.h>

namespace mi
{
namespace opengl
{

Viewport::Viewport()
{
}

void Viewport::set(int x, int y, int width, int height)
{
    this->x = x;
    this->y = y;
    if (height == 0)
    {
        height = 1;
    }
    this->width = width;
    this->height = height;
}

int Viewport::getX()
{
    return x;
}

int Viewport::getY()
{
    return y;
}

int Viewport::getWidth()
{
    return width;
}

int Viewport::getHeight()
{
    return height;
}

int*Viewport::toArray()
{
    int *array = new int[4];
    array[0] = x;
    array[1] = y;
    array[2] = width;
    array[3] = height;
    return array;
}

void Viewport::render()
{
    glViewport(x, y, width, height);
}

}
}
