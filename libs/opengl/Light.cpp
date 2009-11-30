#include <opengl/Light.h>

#include <opengl/OpenGL.h>

using namespace mi::math;

namespace mi
{
namespace opengl
{

Light::Light(int lightId, float *ambientColor)
{
    this->lightId = lightId;
    this->ambientColor = new float[4];
    if (ambientColor == NULL)
    {
        this->ambientColor[0] = 0.2f;
        this->ambientColor[1] = 0.2f;
        this->ambientColor[2] = 0.2f;
        this->ambientColor[3] = 1.0f;
    }
    else
    {
        this->ambientColor[0] = ambientColor[0];
        this->ambientColor[1] = ambientColor[1];
        this->ambientColor[2] = ambientColor[2];
        this->ambientColor[3] = ambientColor[3];
    }
    diffuseColor = new float[4];
    diffuseColor[0] = 1.0f;
    diffuseColor[1] = 1.0f;
    diffuseColor[2] = 1.0f;
    diffuseColor[3] = 1.0f;
    specularColor = new float[4];
    specularColor[0] = 1.0f;
    specularColor[1] = 1.0f;
    specularColor[2] = 1.0f;
    specularColor[3] = 1.0f;
    position = new float[4];
    position[0] = 0.0f;
    position[1] = 0.0f;
    position[2] = 0.0f;
    position[3] = 1.0f;
}

void Light::setPosition(Vector4<float> &position)
{
    this->position[0] = position.getX();
    this->position[1] = position.getY();
    this->position[2] = position.getZ();
    this->position[3] = position.getW();
}

void Light::setAmbientColor(float color0, float color1, float color2, float color3)
{
    ambientColor[0] = color0;
    ambientColor[1] = color1;
    ambientColor[2] = color2;
    ambientColor[3] = color3;
}

void Light::setDiffuseColor(float color0, float color1, float color2, float color3)
{
    diffuseColor[0] = color0;
    diffuseColor[1] = color1;
    diffuseColor[2] = color2;
    diffuseColor[3] = color3;
}

void Light::setSpecularColor(float color0, float color1, float color2, float color3)
{
    specularColor[0] = color0;
    specularColor[1] = color1;
    specularColor[2] = color2;
    specularColor[3] = color3;
}

void Light::render()
{
    glLightfv(lightId, GL_AMBIENT, ambientColor);
    glLightfv(lightId, GL_DIFFUSE, diffuseColor);
    glLightfv(lightId, GL_SPECULAR, specularColor);
    glLightfv(lightId, GL_POSITION, position);
}

Vector4<float> Light::getPosition()
{
    return Vector4<float>(position[0], position[1], position[2], position[3]);
}

float*Light::getAmbientColor()
{
    return ambientColor;
}

float*Light::getDiffuseColor()
{
    return diffuseColor;
}

float*Light::getSpecularColor()
{
    return specularColor;
}

void Light::glEnable()
{
    ::glEnable(lightId);
}

void Light::glDisable()
{
    ::glDisable(lightId);
}

}
}
