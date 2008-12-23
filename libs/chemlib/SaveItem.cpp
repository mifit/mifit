#include "SaveItem.h"
#include "model.h"

namespace chemlib {

SaveItem::SaveItem() {
  SaveMolecule = NULL;
  Title = "Error: empty item";
}

SaveItem::SaveItem(MIMoleculeBase* node, std::string title) {
  SaveMolecule = node;
  Title = title;
}

}
