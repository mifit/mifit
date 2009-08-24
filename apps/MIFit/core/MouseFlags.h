#ifndef MACWINDOWS
#define MACWINDOWS

// mouse messages
#ifndef MK_LBUTTON
#define MK_LBUTTON  0x0001
#endif
#ifndef MK_RBUTTON
#define MK_RBUTTON  0x0002
#endif
#ifndef MK_SHIFT
#define MK_SHIFT    0x0004
#endif
#ifndef MK_CONTROL
#define MK_CONTROL  0x0008
#endif
#ifndef MK_MBUTTON
#define MK_MBUTTON  0x0010
#endif
#undef MK_ALT
#define MK_ALT      0x0020
#define MK_META     0x0040

#endif
