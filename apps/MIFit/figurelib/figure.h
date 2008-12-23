#ifndef FIGURE_H
#define FIGURE_H

#include "chemlib.h"


namespace moldraw {

bool GenerateFigureASP(Drawing *dp,
                       chemlib::MIMoleculeBase* mol,
                       const char* lig_resnumber,
                       const unsigned short chain_id);


} //namespace moldraw
#endif //FIGURE_H
