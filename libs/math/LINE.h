#ifndef MIFIT_MODEL_LINE_H_
#define MIFIT_MODEL_LINE_H_

/**
 * A line between two points.
 */
class LINE
{
public:
    float x1;
    float y1;
    float z1;
    float x2;
    float y2;
    float z2;
    short color;
    short color2;
    /**
     * The atom the line came from
     */
    void *from;
};

#endif /*MIFIT_MODEL_LINE_H_*/
