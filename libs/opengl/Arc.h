#ifndef mi_opengl_Arc_h
#define mi_opengl_Arc_h

#include <math/Vector3.h>

namespace mi
{
    namespace opengl
    {

        class Arc
        {

            float *vertexBuffer;

            int vertexIndex;

            int divisions;

            int bufferSize;

            mi::math::Vector3<float> from;

            mi::math::Vector3<float> to;

            void bisect(mi::math::Vector3<float> &v0, mi::math::Vector3<float> &v1, int level);

            void clearVertices();

            void createVertices();

        public:

            Arc();

            Arc(mi::math::Vector3<float> &from, mi::math::Vector3<float> &to, int divisions);

            mi::math::Vector3<float> getFrom();

            void setFrom(mi::math::Vector3<float> &from);

            mi::math::Vector3<float> getTo();

            void setTo(mi::math::Vector3<float> &to);

            int getDivisions();

            void setDivisions(int divisions);

            void render();

        };

    }
}

#endif // ifndef mi_opengl_Arc_h
