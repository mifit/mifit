#ifndef mi_opengl_interact_TranslationFeedback_h
#define mi_opengl_interact_TranslationFeedback_h

#include <math/Vector3.h>
#include <math/Point3.h>

namespace mi
{
    namespace opengl
    {
        namespace interact
        {

            class TranslationFeedback
            {

                mi::math::Vector3<float> beginPosition;

                mi::math::Vector3<float> endPosition;

                mi::math::Point3<float> color;

                float length;

            public:

                TranslationFeedback(float *color);

                void setLength(float length);

                void render();

                void setFrom(const mi::math::Vector3<float> &from);

                void setTo(const mi::math::Vector3<float> &to);

            };

        }
    }
}

#endif // ifndef mi_opengl_interact_TranslationFeedback_h
