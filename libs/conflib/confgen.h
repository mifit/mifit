#include "chemlib.h"

namespace conflib {


void GenerateDictionary(chemlib::RESIDUE* res,
                        std::vector<chemlib::Bond>& bonds,
                        std::vector<chemlib::Bond>& BondLengths,
                        std::vector<chemlib::ANGLE>& Angles,
                        std::vector<chemlib::TORSION>& Torsions,
                        std::vector<chemlib::TORSION>& Impropers,
                        std::vector<chemlib::PLANE>& Planes,
                        std::vector<chemlib::CHIRAL>& Chirals);

int GenerateEnsemble(chemlib::RESIDUE* res,
                     std::vector<chemlib::Bond>& bonds,
                     chemlib::MIMolDictionary* dictionary,
                     bool replace = false);

int GenerateEnsemble(chemlib::RESIDUE* res,
                     chemlib::MIMoleculeBase* model,
                     std::vector<chemlib::Bond>& bonds,
                     std::vector<chemlib::TORSION>& torsions,
                     chemlib::GeomSaver& confs);

void GenerateCoordinates(chemlib::RESIDUE* res,
                         const std::vector<chemlib::Bond>& bonds,
                         std::string& log);

void GenerateCoordinates(chemlib::Ligand* lig, std::string& log);

void sdgGenerateCoordinates(chemlib::Ligand* lig);

void sdgRefine(chemlib::Ligand* lig);
}
