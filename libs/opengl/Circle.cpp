#include <opengl/Circle.h>

#include <opengl/OpenGL.h>
#include <math/mathlib.h>
#include <cmath>
#include <cstdlib>

namespace mi
{
namespace opengl
{

Circle::Circle(float radius)
{
    setRadius(radius);
    setAngle((float) (2.0 * M_PI));
    setDetail(80);
}

Circle::Circle(float radius, float angle)
{
    setRadius(radius);
    setAngle(angle);
    setDetail(80);
}

float Circle::getRadius()
{
    return radius;
}

void Circle::setRadius(float radius)
{
    radius = std::abs(radius);
    if (this->radius != radius)
    {
        this->radius = radius;
    }
}

int Circle::getDetail()
{
    return detail;
}

void Circle::setDetail(int detail)
{
    detail = std::abs(detail);
    if (detail < 6)
    {
        detail = 6;
    }
    if (detail != this->detail)
    {
        this->detail = detail;
    }
}

float Circle::getAngle()
{
    return angle;
}

void Circle::setAngle(float angle)
{
    if (this->angle != angle)
    {
        this->angle = angle;
    }
}

void Circle::render()
{

    glPushAttrib(GL_ENABLE_BIT);
    glDisable(GL_TEXTURE_2D);
    glDisable(GL_LIGHTING);
    glDisable(GL_BLEND);

    glBegin(GL_LINE_STRIP);

    for (int i = 0; i <= detail; i++)
    {
        float x = radius * (float) std::cos(i * angle / detail);
        float y = radius * (float) std::sin(i * angle / detail);
        glVertex3f(x, y, 0.0f);
    }

    glEnd();
    glPopAttrib();
}

}
}
