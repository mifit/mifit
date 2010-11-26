#include "ViewPointIO.h"

#include "Cfiles.h"
#include "ViewPoint.h"

void ViewPointIO::save(ViewPoint &vp, CArchive &ar)
{
    if (ar.IsStoring())
    {
        std::string s;
        s = format("translation %0.2f %0.2f %0.2f\n",
                   vp.getcenter(0), vp.getcenter(1), vp.getcenter(2));
        ar.Write(s.c_str(), s.size());
        float viewmat[3][3];
        vp.copymatrix(viewmat);
        s = format("rotation %0.4f %0.4f %0.4f %0.4f %0.4f %0.4f %0.4f %0.4f %0.4f\n",
                   viewmat[0][0], viewmat[0][1], viewmat[0][2],
                   viewmat[1][0], viewmat[1][1], viewmat[1][2],
                   viewmat[2][0], viewmat[2][1], viewmat[2][2]);
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
                float viewmat[3][3];
                viewmat[0][0] = v1;
                viewmat[0][1] = v2;
                viewmat[0][2] = v3;
                viewmat[1][0] = v4;
                viewmat[1][1] = v5;
                viewmat[1][2] = v6;
                viewmat[2][0] = v7;
                viewmat[2][1] = v8;
                viewmat[2][2] = v9;
                vp.setmatrix(viewmat);
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
                vp.moveto(v1, v2, v3);
            }
        }
    }
}

