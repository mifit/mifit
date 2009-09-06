#ifndef mifit_ui_JobsView_h
#define mifit_ui_JobsView_h

#include <QWidget>

class BatchJob;
class BatchJobManager;
class JobsTree;

class QTreeWidgetItem;

/**
 * Class to control the user interface view of jobs.
 */
class JobsView : public QWidget {

  JobsTree* jobsTree;
  bool listeningToBatchJobManager;

  void addJobToTree(BatchJob* job);
  bool jobInTree(BatchJob* job);
  void stylizeItem(QTreeWidgetItem* id, BatchJob* job);

public slots:
  void jobAdded(BatchJob* job);
  void jobDeleted(BatchJob* job);
  void jobChanged(BatchJob* job);

public:

  JobsView(QWidget* parent);
  ~JobsView();

  void update(BatchJobManager* batchJobManager);
};

#endif
