#include <vector>
#include "UnifiedRes.h"
#include "Shape.h"
#include "Sun.h"


using namespace moldraw;

void UnifiedRes::CalcDirection(double* direction, double* extent) {
  if (partners.empty()) {
    *direction = 0;
    *extent = 360.0;
  } else {

    *direction = atan2(partners.front()->x() - x,
                   partners.front()->y() - y);
    *extent = 120.0;
  }
}

void UnifiedRes::Draw(Drawing* dp) {
  double direction, extent;
  CalcDirection(&direction, &extent);

  std::string label = type;
  label += " ";
  label += name;

  if (chain_id != ' ') {
    label += "(";
    label += chain_id;
    label += ")";
  }

  Sun s(dp,
        x,
        y,
        Drawing::UNIRES_RADIUS,
        sin(direction),
        cos(direction),
        extent,
        Drawing::SUN_LABEL_SIZE,
        Drawing::SUN_BORDER_WIDTH,
        label);

  s.Draw();
}

