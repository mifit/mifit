#include "ViewPointIO.h"

#include "Cfiles.h"
#include "ViewPoint.h"
#include <math/Matrix4.h>

using namespace mi::math;

void ViewPointIO::save(ViewPoint &vp, CArchive &ar)
{
    if (ar.IsStoring())
    {
        std::string s;
        s = format("translation %0.2f %0.2f %0.2f\n",
                   vp.center()[0], vp.center()[0], vp.center()[2]);
        ar.Write(s.c_str(), s.size());
        Matrix4<float> viewmat(vp.orientation());
        s = format("rotation %0.4f %0.4f %0.4f %0.4f %0.4f %0.4f %0.4f %0.4f %0.4f\n",
                   viewmat.m00, viewmat.m01, viewmat.m02,
                   viewmat.m10, viewmat.m11, viewmat.m12,
                   viewmat.m20, viewmat.m21, viewmat.m22);
        ar.Write(s.c_str(), s.size());
        s = format("zoom %0.3f\nperspective %0.3f\n", vp.scale(), vp.perspective());
        ar.Write(s.c_str(), s.size());
        s = format("frontclip %0.2f\nbackclip %0.2f\n", vp.frontClip(), vp.backClip());
        ar.Write(s.c_str(), s.size());
        s = format("transform\n");
        ar.Write(s.c_str(), s.size());
    }
}

void ViewPointIO::load(ViewPoint &vp, FILE *fp)
{
    float v1, v2, v3, v4, v5, v6, v7, v8, v9;
    std::string buf;

    std::auto_ptr<io> ioObj(io::defaultIo());
    io &file = *ioObj;
    file.attach(fp);
    while (file.readLine(buf) != 0)
    {
        std::transform(buf.begin(), buf.end(), buf.begin(), tolower);
        if (!strncmp(buf.c_str(), "rotation", 8))
        {
            if (sscanf(buf.c_str(), "%*s%f%f%f%f%f%f%f%f%f", &v1, &v2, &v3, &v4, &v5, &v6, &v7, &v8, &v9) == 9)
            {
                vp.setView(Matrix4<float>(v1, v2, v3, 0.0f,
                                      v4, v5, v6, 0.0f,
                                      v7, v8, v9, 0.0f,
                                      0.0f, 0.0f, 0.0f, 1.0f));
            }
        }
        else if (!strncmp(buf.c_str(), "zoom", 4))
        {
            if (sscanf(buf.c_str(), "%*s%f", &v1) == 1)
            {
                vp.setScale(ROUND(v1));
            }
        }
        else if (!strncmp(buf.c_str(), "frontclip", 9))
        {
            if (sscanf(buf.c_str(), "%*s%f", &v1) == 1)
            {
                vp.setFrontClip(v1);
            }
        }
        else if (!strncmp(buf.c_str(), "backclip", 8))
        {
            if (sscanf(buf.c_str(), "%*s%f", &v1) == 1)
            {
                vp.setBackClip(v1);
            }
        }
        else if (!strncmp(buf.c_str(), "perspect", 8))
        {
            if (sscanf(buf.c_str(), "%*s%f", &v1) == 1)
            {
                vp.setPerspective(v1);
            }
        }
        else if (!strncmp(buf.c_str(), "translation", 11))
        {
            if (sscanf(buf.c_str(), "%*s%f%f%f", &v1, &v2, &v3) == 3)
            {
                vp.setCenter(mi::math::Vector3<float>(v1, v2, v3));
            }
        }
    }
}

