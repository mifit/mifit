#include "JobsView.h"

#include <QApplication>
#include <QVBoxLayout>

#include <vector>
#include <map>

#include <nongui/nonguilib.h>
#include "jobs/jobslib.h"
#include "core/MIData.h"

#include "MIEventHandler.h"
#include "MIEventHandlerMacros.h"
#include "MIMainWindow.h"
#include "MIMenu.h"
#include "MIQTreeWidget.h"

#include "uitest.h"
#include "ui/MIDialog.h"
#include "id.h"

#include <images/jobsList.xpm>
#include <images/job.xpm>
#include <images/jobRunning.xpm>
#include <images/jobOk.xpm>
#include <images/jobError.xpm>




class JobsTree : public MIQTreeWidget, public MIEventHandler {
Q_OBJECT

public:
  JobsTree(QWidget* parent);
  virtual ~JobsTree();

 private slots:
  void OnItemClicked(QTreeWidgetItem *item, int column); // single click
  void OnItemActivated(QTreeWidgetItem *item, int column); // double click
  void OnItemPressed(QTreeWidgetItem *item, int column); // possible right click

  void DeleteJob();
  void ShowProperties();
  void ShowLog();
  void CleanSuccessful();
  void CleanAll();
  void DetachJob();
  void OpenResults();

public:
  QTreeWidgetItem* rootId;
  BatchJobManager* batchJobManager;
  std::map<BatchJob*, QTreeWidgetItem*> jobToItem;
  std::map<QTreeWidgetItem*, BatchJob*> itemToJob;
};



JobsTree::JobsTree(QWidget* parent): MIQTreeWidget(parent), MIEventHandler(this) {

  std::vector<QIcon> imageList;
  QIcon jobsListImage=QIcon(QPixmap(jobsList_xpm));
  imageList.push_back(jobsListImage);
  QIcon jobImage(job_xpm);
  imageList.push_back(jobImage);
  QIcon jobRunningImage(jobRunning_xpm);
  imageList.push_back(jobRunningImage);
  QIcon jobOkImage(jobOk_xpm);
  imageList.push_back(jobOkImage);
  QIcon jobErrorImage(jobError_xpm);
  imageList.push_back(jobErrorImage);
  AssignImageList(imageList);

  std::string rootText = std::string("Jobs List");
  setHeaderLabel(rootText.c_str());
  rootId = invisibleRootItem();

  connect(this, SIGNAL(itemPressed(QTreeWidgetItem *, int)),
          this, SLOT(OnItemPressed(QTreeWidgetItem *, int)));

  connect(this, SIGNAL(itemClicked(QTreeWidgetItem *, int)),
          this, SLOT(OnItemClicked(QTreeWidgetItem *, int)));

  connect(this, SIGNAL(itemActivated(QTreeWidgetItem *, int)),
          this, SLOT(OnItemActivated(QTreeWidgetItem *, int)));

BEGIN_EVENT_TABLE(this, none)
EVT_MENU(ID_JOBSVIEW_DELETE, JobsTree::DeleteJob)
EVT_MENU(ID_JOBSVIEW_PROPERTIES, JobsTree::ShowProperties)
EVT_MENU(ID_JOBSVIEW_SHOWLOG, JobsTree::ShowLog)
EVT_MENU(ID_JOBSVIEW_OPENRESULTS, JobsTree::OpenResults)
EVT_MENU(ID_JOBSVIEW_CLEANSUCCESSFUL, JobsTree::CleanSuccessful)
EVT_MENU(ID_JOBSVIEW_CLEANALL, JobsTree::CleanAll)
EVT_MENU(ID_JOBSVIEW_DETACH, JobsTree::DetachJob)
END_EVENT_TABLE()
}

JobsTree::~JobsTree() {
}

void JobsTree::OnItemClicked(QTreeWidgetItem *, int) {
}

void JobsTree::OnItemActivated(QTreeWidgetItem *, int) {
  QList<QTreeWidgetItem *> selected;
  GetSelections(selected);
  if (selected.size() == 0) {
    return;
  }
  for (int i = 0; i < selected.size(); i++) {
    QTreeWidgetItem* id = selected[i];
    BatchJob* job = itemToJob[id];
    job->openResults();
  }
}

void JobsTree::OnItemPressed(QTreeWidgetItem *id, int) {
  if (QApplication::mouseButtons() != Qt::RightButton)
    return;

  QList<QTreeWidgetItem *> selected;
  GetSelections(selected);
  if (selected.size() <= 1) {
    //set selected item to right-clicked item
    clearSelection();
    id->setSelected(true);
  }
  MIMenu* menu = new MIMenu(*this);
  menu->Append(ID_JOBSVIEW_OPENRESULTS, "Open Results...");
  menu->Append(ID_JOBSVIEW_SHOWLOG, "Show Log File", "Open the log file in a window", false);
  menu->Append(ID_JOBSVIEW_PROPERTIES, "Job Properties", "Properties for this job", false);
  menu->Append(ID_JOBSVIEW_DELETE, "Delete Job", "Delete this job", false);
  menu->Append(ID_JOBSVIEW_CLEANSUCCESSFUL, "Clean Successful", "Remove jobs that have completed successfully", false);
  menu->Append(ID_JOBSVIEW_CLEANALL, "Clean All", "Remove jobs that have completed", false);
  menu->Append(ID_JOBSVIEW_DETACH, "Detach Job", "Detach this job", false);
  BatchJob* job = NULL;
  if (itemToJob.find(id) != itemToJob.end()) {
    job = itemToJob[id];
  }
  if (id == rootId) {
    menu->Enable(ID_JOBSVIEW_DELETE, false);
    menu->Enable(ID_JOBSVIEW_PROPERTIES, false);
    menu->Enable(ID_JOBSVIEW_SHOWLOG, false);
    menu->Enable(ID_JOBSVIEW_DETACH, false);
  } else {
    if (job != NULL && job->IsRunning()) {
      menu->Enable(ID_JOBSVIEW_OPENRESULTS, false);
      menu->Enable(ID_JOBSVIEW_SHOWLOG, false);
      menu->Enable(ID_JOBSVIEW_DETACH, true);
    } else {
      menu->Enable(ID_JOBSVIEW_DETACH, false);
    }
  }

  QPoint pos=QCursor::pos();
  menu->doExec(pos);
}

void JobsTree::DeleteJob() {
  QList<QTreeWidgetItem *> selected;
  GetSelections(selected);
  if (selected.size() == 0) {
    return;
  } else if (selected.size() > 1) {
    if (MIMessageBox("Are you sure you want to delete the selected jobs?",
          "Confirm Delete Jobs", MIDIALOG_YES_NO) != MI_YES) {
      return;
    }
  } else {
    QTreeWidgetItem* id = selected[0];
    BatchJob* job = itemToJob[id];
    if (!job) {
      return;
    }
    std::string message;
    message=::format("Are you sure you want to delete job %d?", job->JobId);
    if (MIMessageBox(message.c_str(), "Confirm Delete Job", MIDIALOG_YES_NO) != MI_YES) {
      return;
    }
  }
  for (int i = 0; i < selected.size(); i++) {
    QTreeWidgetItem* id = selected[i];
    BatchJob* job = itemToJob[id];
    batchJobManager->DeleteJob(job);
  }
}

void JobsTree::DetachJob() {
  QList<QTreeWidgetItem *> selected;
  GetSelections(selected);
  if (selected.size() == 0) {
    return;
  } else if (selected.size() > 1) {
    if (MIMessageBox("Are you sure you want to detach the selected jobs?",
          "Confirm Detach Jobs", MIDIALOG_YES_NO) != MI_YES) {
      return;
    }
  } else {
    QTreeWidgetItem* id = selected[0];
    BatchJob* job = itemToJob[id];
    std::string message;
    message=::format("Are you sure you want to detach job %d?", job->JobId);
    if (MIMessageBox(message.c_str(), "Confirm Detach Job", MIDIALOG_YES_NO) != MI_YES) {
      return;
    }
  }
  for (int i = 0; i < selected.size(); i++) {
    QTreeWidgetItem* id = selected[i];
    BatchJob* job = itemToJob[id];
    batchJobManager->DetachJob(job);
  }
}

void JobsTree::ShowProperties() {
  QList<QTreeWidgetItem*> selected;
  GetSelections(selected);
  for (int i = 0; i < selected.size(); i++) {
    QTreeWidgetItem* id = selected[i];
    BatchJob* job = itemToJob[id];
    MIMessageBox(job->Info().c_str(), "Job Properties", MIDIALOG_ICON_INFORMATION);
  }
}

void JobsTree::ShowLog() {
  QList<QTreeWidgetItem*> selected;
  GetSelections(selected);
  for (int i = 0; i < selected.size(); i++) {
    QTreeWidgetItem* id = selected[i];
    BatchJob* job = itemToJob[id];
    batchJobManager->ShowLogFile(job);
  }
}

void JobsTree::CleanSuccessful() {
  if (MIMessageBox("Are you sure you want to delete all successful jobs?",
        "Confirm Delete Jobs", MIDIALOG_YES_NO) != MI_YES) {
    return;
  }
  batchJobManager->CleanSucc();
}

void JobsTree::CleanAll() {
  if (MIMessageBox("Are you sure you want to delete all jobs?",
        "Confirm Delete Jobs", MIDIALOG_YES_NO) != MI_YES) {
    return;
  }
  batchJobManager->CleanAll();
}

void JobsTree::OpenResults() {
  QList<QTreeWidgetItem*> selected;
  GetSelections(selected);
  for (int i = 0; i < selected.size(); i++) {
    QTreeWidgetItem* id = selected[i];
    BatchJob* job = itemToJob[id];
    job->openResults();
  }
}

JobsView::JobsView(QWidget* parent) : QWidget(parent), listeningToBatchJobManager(false) {

  jobsTree = new JobsTree(this);
  QVBoxLayout *vbox=new QVBoxLayout();
  vbox->setContentsMargins(0, 0, 0, 0);
  vbox->setSpacing(2);
  vbox->addWidget(jobsTree);
  setLayout(vbox);
}

JobsView::~JobsView() {
}

void JobsView::jobAdded(BatchJob* job) {
  addJobToTree(job);
  connect(job, SIGNAL(jobChanged(BatchJob*)),
          this, SLOT(jobChanged(BatchJob*)));
}

void JobsView::jobDeleted(BatchJob* job) {
  QTreeWidgetItem* id = jobsTree->jobToItem[job];
  jobsTree->Delete(id);
  jobsTree->jobToItem.erase(job);
  jobsTree->itemToJob.erase(id);
}

void JobsView::jobChanged(BatchJob* job) {
  QTreeWidgetItem* id = jobsTree->jobToItem[job];
  stylizeItem(id, job);
}

void JobsView::update(BatchJobManager* batchJobManager) {
  if (!listeningToBatchJobManager) {
    listeningToBatchJobManager = true;
    connect(batchJobManager, SIGNAL(jobAdded(BatchJob*)),
            this, SLOT(jobAdded(BatchJob*)));
    connect(batchJobManager, SIGNAL(jobDeleted(BatchJob*)),
            this, SLOT(jobDeleted(BatchJob*)));
    jobsTree->batchJobManager = batchJobManager;
  }
  std::vector<BatchJob*>& jobList = *batchJobManager->GetJobList();
  std::vector<BatchJob*>::iterator jobIter = jobList.begin();
  while (jobIter != jobList.end()) {
    BatchJob* job = *jobIter;
    ++jobIter;
    if (!jobInTree(job)) {
      addJobToTree(job);
    }
  }
}

bool JobsView::jobInTree(BatchJob* job) {
  return jobsTree->jobToItem.find(job) != jobsTree->jobToItem.end();
}

void JobsView::addJobToTree(BatchJob* job) {
  QTreeWidgetItem* id = jobsTree->appendItem(jobsTree->rootId, "New Job", -1, -1);
  jobsTree->jobToItem[job] = id;
  jobsTree->itemToJob[id] = job;
  stylizeItem(id, job);
}

void JobsView::stylizeItem(QTreeWidgetItem* id, BatchJob* job) {
  int image = 1;
  if (job->IsRunning()) {
    image = 2;
  } else if (job->isCompleted()) {
    if (job->RanOK()) {
      image = 3;
    } else {
      image = 4;
    }
  }
  id->setIcon(0,jobsTree->GetIcon(image));
  std::string jobText;
  if (job->getSettings()["jobName"].str != MIDatum::INVALID_STRING) {
    jobText=::format("%s (%d)", job->getSettings()["jobName"].str.c_str(), job->JobId);
  } else {
    jobText=::format("%d", job->JobId);
  }
  id->setText(0, jobText.c_str());
}


#include "JobsView.moc"
