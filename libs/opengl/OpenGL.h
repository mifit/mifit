#ifndef mi_opengl_OpenGL_h
#define mi_opengl_OpenGL_h

#ifdef _WIN32
#include <windows.h>
#undef min
#undef max
#endif

#if !defined(__APPLE__) && !defined(__WXCOCOA__)
#include <GL/gl.h>
#include <GL/glu.h>
#ifdef _WIN32
#include <GL/glext.h>
#endif
#else
#include <OpenGL/gl.h>
#include <OpenGL/glu.h>
#endif

#endif
