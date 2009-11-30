#ifndef MIFIT_MODEL_CREFL_H_
#define MIFIT_MODEL_CREFL_H_

/**
 * A refection structure.
 */
class CREFL
{
public:

    int ind[3];

    /**
     * as read in
     */
    float fo;

    /**
     * as read in
     */
    float sigma;

    /**
     * as read in
     */
    float fc;

    /**
     * as read in
     */
    float phi;

    float sthol;

    /**
     * coefficients of current map as set in nextf in fftsub.c
     */
    float coef;

    float fom;

    /**
     * a,b of current calculation
     */
    float acalc;

    /**
     * a,b of current calculation
     */
    float bcalc;

    /**
     * a,b of all atoms
     */
    float awhole;

    /**
     * a,b of all atoms
     */
    float bwhole;

    /**
     * a,b of non-current structure
     */
    float astatic;

    /**
     * a,b of non-current structure
     */
    float bstatic;

    /**
     * a,b of current part of structure
     */
    float apart;

    /**
     * a,b of current part of structure
     */
    float bpart;

    /**
     * used to store R-free flag info
     */
    short freeRflag;
};


#endif /*MIFIT_MODEL_CREFL_H_*/
