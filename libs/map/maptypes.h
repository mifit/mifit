#ifndef MI_MAPTYPES_H
#define MI_MAPTYPES_H

namespace MIMapType
{
    const unsigned int Fo = 0;
    const unsigned int Fc = 1;
    const unsigned int TwoFoFc = 2;
    const unsigned int FoFc = 3;
    const unsigned int FoFo = 4;
    const unsigned int Fofom = 5;
    const unsigned int ThreeFoTwoFc = 6;
    const unsigned int FiveFoThreeFc = 7;
    const unsigned int TwoFoFcSigmaA = 8;
    const unsigned int FoFcSigmaA = 9;
    const unsigned int DirectFFT = 10; // has same effect as Fo!
}

/*  prime is largest prime number accepted by fft, even_odd is
 *  whether number must be even or not.
 *    1 is odd or even , 2 is even only
 *    note when using fft3d nx must be even, ny and nz can be either
 *    odd or even -  prime should be 7 */
#define FFT_PRIME 7
#define EVEN 2
#define ODDOREVEN 1

const char *StringForMapType(unsigned int maptype);
unsigned int MapTypeForString(const std::string &str);

// returns true if this is generally a type for a regular map, false if diff map
bool IsRegularMapType(unsigned int maptype);
#endif // ifndef MI_MAPTYPES_H
