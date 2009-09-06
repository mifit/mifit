#include <opengl/Sphere.h>

#include <opengl/OpenGL.h>
#include <cmath>
#include <math/mi_math.h>

using namespace std;

namespace mi {
namespace opengl {

Sphere::Sphere(float radius, int slices, int stacks) {

  if (radius <= 0.0f) {
    throw "sphere radius must be > 0.0";
  }
  if (slices <= 2) {
    throw "sphere slices must be > 2";
  }
  if (stacks <= 1) {
    throw "sphere stacks must be > 1";
  }

  vertexCount = 2 * stacks * (slices + 1);
  vertexBuffer = new float[3 * vertexCount];
  normalBuffer = new float[3 * vertexCount];
  textureBuffer = new float[2 * vertexCount];

  float drho = (float) (M_PI) / (float) stacks;
  float deltaTheta = 2.0f * (float) (M_PI) / (float) slices;
  float deltaSlice = 1.0f / (float) slices;
  float deltaStack = 1.0f / (float) stacks;
  float t = 1.0f;

  int vertexIndex = 0;
  int normalIndex = 0;
  int textureIndex = 0;
  for (int i = 0; i < stacks; i++) {
    float rho = (float) i * drho;
    float sinRho = (float) (std::sin(rho));
    float cosRho = (float) (std::cos(rho));
    float srhodrho = (float) (std::sin(rho + drho));
    float crhodrho = (float) (std::cos(rho + drho));

    float s = 0.0f;
    for (int j = 0; j <= slices; j++) {
      float theta = (j == slices) ? 0.0f : j * deltaTheta;
      float sinTheta = (float) -std::sin(theta);
      float cosTheta = (float) std::cos(theta);

      float x = sinTheta * sinRho;
      float y = cosTheta * sinRho;
      float z = cosRho;

      textureBuffer[textureIndex] = s;
      ++textureIndex;
      textureBuffer[textureIndex] = t;
      ++textureIndex;
      normalBuffer[normalIndex] = x;
      ++normalIndex;
      normalBuffer[normalIndex] = y;
      ++normalIndex;
      normalBuffer[normalIndex] = z;
      ++normalIndex;
      vertexBuffer[vertexIndex] = x * radius;
      ++vertexIndex;
      vertexBuffer[vertexIndex] = y * radius;
      ++vertexIndex;
      vertexBuffer[vertexIndex] = z * radius;
      ++vertexIndex;

      x = sinTheta * srhodrho;
      y = cosTheta * srhodrho;
      z = crhodrho;
      textureBuffer[textureIndex] = s;
      ++textureIndex;
      textureBuffer[textureIndex] = t - deltaStack;
      ++textureIndex;
      s += deltaSlice;
      normalBuffer[normalIndex] = x;
      ++normalIndex;
      normalBuffer[normalIndex] = y;
      ++normalIndex;
      normalBuffer[normalIndex] = z;
      ++normalIndex;
      vertexBuffer[vertexIndex] = x * radius;
      ++vertexIndex;
      vertexBuffer[vertexIndex] = y * radius;
      ++vertexIndex;
      vertexBuffer[vertexIndex] = z * radius;
      ++vertexIndex;

    }

    t -= deltaStack;
  }

}

void Sphere::render() {
  glEnableClientState(GL_VERTEX_ARRAY);
  glEnableClientState(GL_NORMAL_ARRAY);
  glEnableClientState(GL_TEXTURE_COORD_ARRAY);

  glVertexPointer(3, GL_FLOAT, 0, vertexBuffer);
  glNormalPointer(GL_FLOAT, 0, normalBuffer);
  glTexCoordPointer(2, GL_FLOAT, 0, textureBuffer);

  glDrawArrays(GL_TRIANGLE_STRIP, 0, vertexCount);
}

}
}
