#ifndef MIFIT_UI_COLORUTIL_H_
#define MIFIT_UI_COLORUTIL_H_

short secstrcolor(char);
short bvaluecolor(float);

//@{
// Init the internal color names.
// Run upon program start.
//@}
void init_colornames();

//@{
// returns a color (see colors.h) given an atom name (such as CB).
//@}
short color_by_name(const char*);
//@{
// returns a color (see colors.h) given an atom name (such as CB) and an altloc character.
//@}
short color_by_name(const char*, char);

#endif /*MIFIT_UI_COLORUTIL_H_*/
