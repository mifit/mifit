#ifndef mifit_map_MINATOM_h
#define mifit_map_MINATOM_h

/**
 * A struct for saving the location of an atom with minimal overhead.
 * Used in the find water function.
 */
typedef struct minatom
{
    float fx, fy, fz;
    float tx, ty, tz;
    float cx, cy, cz;
    char type;
    int symm;
} MINATOM;

#endif
