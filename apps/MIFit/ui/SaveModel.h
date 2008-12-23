#ifndef MIFIT_UI_SAVEMODEL_H_
#define MIFIT_UI_SAVEMODEL_H_

/**
 * Save model information for the checkpoint operation.
 */
class SaveModel {
public:
  /**
   * The file the model has been checkpointed into.
   */
  std::string path;

  /**
   * The model is from.
   */
  Molecule* model;

  /**
   * The date and time of the checkpointing.
   */
  std::string time;

  /**
   * An informative description so the user can figure out where the best checkpointed model is.
   */
  std::string description;
};

#endif /*MIFIT_UI_SAVEMODEL_H_*/
