#include "corelib.h"

#include "SurfaceSphere.h"


SurfaceSphere::SurfaceSphere() :  spacing(-1.0f), radius(-1.0f) {
}

float SurfaceSphere::getRadius() {
  return radius;
}

float SurfaceSphere::getSpacing() {
  return spacing;
}

std::vector<APOINT>& SurfaceSphere::getPoints() {
  return points;
}

void SurfaceSphere::clearPoints() {
  std::vector<APOINT>().swap(points); // was points.clear();
  radius = -1.0f;
  spacing = -1.0f;
}

void SurfaceSphere::build(float rad, float spac) {
  if (spac < 0.0f) {
    return;
  }
  if (rad <= 0.0f) {
    return;
  }
  Logger::log("Creating sphere of radius %f spacing %f\n",rad,spac);
  spacing = spac;
  radius = rad;
  int color = Colors::WHITE;
  float pi = (float) M_PI;
  float twopi = (float) M_PI * 2.0F;

  float a = pi*radius/spacing +0.5F;
  int nlevels = (int)a;

  points.clear();
  for (int i = 0; i < nlevels; i++) {
    double psi = pi * (float)i/(float)nlevels;
    float rcirc = radius * (float)std::sin(psi);
    float z = radius * (float)std::cos(psi);
    a =  twopi*rcirc/spacing +0.5F;
    int ncirc = (int)a;
    for (int j = 0; j < ncirc; j++) {
      APOINT p;
      a = (float)j/(float)ncirc;
      float phi = twopi*a;
      p.z = z;
      p.x = rcirc*(float)std::cos(phi);
      p.y = rcirc*(float)std::sin(phi);
      p.color = color;
      points.push_back(p);
    }
  }
  // Add the poles
  APOINT pole1;
  pole1.z = -radius;
  pole1.x = 0.0;
  pole1.y = 0.0;
  pole1.color = color;
  points.push_back(pole1);
  APOINT pole2;
  pole2.z = radius;
  pole2.x = 0.0;
  pole2.y = 0.0;
  pole2.color = color;
  points.push_back(pole2);
}

