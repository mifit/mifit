#include <opengl/interact/TranslationFeedback.h>

#include <opengl/OpenGL.h>

using namespace mi::math;
using namespace mi::opengl;

namespace mi
{
namespace opengl
{
namespace interact
{

TranslationFeedback::TranslationFeedback(float *color)
{
    this->color.set(color[0], color[1], color[2]);
}

void TranslationFeedback::setLength(float length)
{
    this->length = length;
}

void TranslationFeedback::render()
{

    glPushAttrib(GL_DEPTH_BUFFER_BIT | GL_CURRENT_BIT | GL_LIGHTING_BIT | GL_LINE_BIT);
    glEnable(GL_DEPTH_TEST);
    glDisable(GL_LIGHTING);
    glLineWidth(1.0f);

    glColor3f(color.getX(), color.getY(), color.getZ());

    glBegin(GL_LINES);
    glVertex3f(beginPosition.getX(), beginPosition.getY(), beginPosition.getZ());
    glVertex3f(endPosition.getX(), endPosition.getY(), endPosition.getZ());
    glEnd();

    glPushMatrix();
    glTranslatef(beginPosition.getX(), beginPosition.getY(), beginPosition.getZ());

    glBegin(GL_LINES);

    float margin = length / 3;
    glVertex3d(-length, 0, 0);
    glVertex3d(-margin, 0, 0);
    glVertex3d(margin, 0, 0);
    glVertex3d(length, 0, 0);

    glVertex3d(0, -length, 0);
    glVertex3d(0, -margin, 0);
    glVertex3d(0, margin, 0);
    glVertex3d(0, length, 0);

    glVertex3d(0, 0, -length);
    glVertex3d(0, 0, -margin);
    glVertex3d(0, 0, margin);
    glVertex3d(0, 0, length);

    glEnd();

    glPopMatrix();

    glPopAttrib();

}

void TranslationFeedback::setFrom(const Vector3<float> &from)
{
    beginPosition.set(from);
}

void TranslationFeedback::setTo(const Vector3<float> &to)
{
    endPosition.set(to);
}

}
}
}
