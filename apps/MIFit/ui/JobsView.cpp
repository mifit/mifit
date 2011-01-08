#include "JobsView.h"

#include <QApplication>
#include <QMenu>
#include <QMessageBox>
#include <QVBoxLayout>

#include <vector>
#include <map>

#include "jobs/jobslib.h"

#include "MIMainWindow.h"
#include "MIQTreeWidget.h"

#include <images/jobsList.xpm>
#include <images/job.xpm>
#include <images/jobRunning.xpm>
#include <images/jobOk.xpm>
#include <images/jobError.xpm>




class JobsTree
    : public MIQTreeWidget
{
    Q_OBJECT

public:
    JobsTree(QWidget *parent);
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
    void OpenResults();

public:
    QTreeWidgetItem *rootId;
    BatchJobManager *batchJobManager;
    std::map<BatchJob*, QTreeWidgetItem*> jobToItem;
    std::map<QTreeWidgetItem*, BatchJob*> itemToJob;
};



JobsTree::JobsTree(QWidget *parent)
    : MIQTreeWidget(parent)
{

    std::vector<QIcon> imageList;
    QIcon jobsListImage = QIcon(QPixmap(jobsList_xpm));
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

}

JobsTree::~JobsTree()
{
}

void JobsTree::OnItemClicked(QTreeWidgetItem*, int)
{
}

void JobsTree::OnItemActivated(QTreeWidgetItem*, int)
{
    QList<QTreeWidgetItem*> selected;
    GetSelections(selected);
    if (selected.size() == 0)
    {
        return;
    }
    for (int i = 0; i < selected.size(); i++)
    {
        QTreeWidgetItem *id = selected[i];
        BatchJob *job = itemToJob[id];
        job->openResults();
    }
}

void JobsTree::OnItemPressed(QTreeWidgetItem *id, int)
{
    if (QApplication::mouseButtons() != Qt::RightButton)
        return;

    QList<QTreeWidgetItem*> selected;
    GetSelections(selected);
    if (selected.size() <= 1)
    {
        //set selected item to right-clicked item
        clearSelection();
        id->setSelected(true);
    }
    QMenu *menu = new QMenu;
    QAction *openResultsAction = menu->addAction("Open Results...", this, SLOT(OpenResults()));
    QAction *showLogAction = menu->addAction("Show Log File", this, SLOT(ShowLog()));
    QAction *propertiesAction = menu->addAction("Job Properties", this, SLOT(ShowProperties()));
    QAction *deleteAction = menu->addAction("Delete Job", this, SLOT(DeleteJob()));
    menu->addAction("Clean Successful", this, SLOT(CleanSuccessful()));
    menu->addAction("Clean All", this, SLOT(CleanAll()));

    BatchJob *job = NULL;
    if (itemToJob.find(id) != itemToJob.end())
        job = itemToJob[id];

    if (id == rootId)
    {
        deleteAction->setEnabled(false);
        propertiesAction->setEnabled(false);
        showLogAction->setEnabled(false);
    }
    else
    {
        if (job != NULL && job->isRunning())
        {
            openResultsAction->setEnabled(false);
            showLogAction->setEnabled(false);
        }
    }

    QPoint pos = QCursor::pos();
    menu->exec(pos);
    delete menu;
}

void JobsTree::DeleteJob()
{
    QList<QTreeWidgetItem*> selected;
    GetSelections(selected);
    if (selected.size() == 0)
    {
        return;
    }
    else if (selected.size() > 1)
    {
        if (QMessageBox::question(0, "Confirm Delete Jobs",
                                  "Are you sure you want to delete the selected jobs?",
                                  QMessageBox::Yes | QMessageBox::No, QMessageBox::No) == QMessageBox::No)
        {
            return;
        }
    }
    else
    {
        QTreeWidgetItem *id = selected[0];
        BatchJob *job = itemToJob[id];
        if (!job)
        {
            return;
        }
        QString message = QString("Are you sure you want to delete job %1?").arg(job->jobId());
        if (QMessageBox::question(0, "Confirm Delete Job",
                                  message,
                                  QMessageBox::Yes | QMessageBox::No, QMessageBox::No) == QMessageBox::No)
        {
            return;
        }
    }
    for (int i = 0; i < selected.size(); i++)
    {
        QTreeWidgetItem *id = selected[i];
        BatchJob *job = itemToJob[id];
        batchJobManager->DeleteJob(job);
    }
}

void JobsTree::ShowProperties()
{
    QList<QTreeWidgetItem*> selected;
    GetSelections(selected);
    for (int i = 0; i < selected.size(); i++)
    {
        QTreeWidgetItem *id = selected[i];
        BatchJob *job = itemToJob[id];
        QMessageBox::information(this, "Job Properties", job->Info());
    }
}

void JobsTree::ShowLog()
{
    QList<QTreeWidgetItem*> selected;
    GetSelections(selected);
    for (int i = 0; i < selected.size(); i++)
    {
        QTreeWidgetItem *id = selected[i];
        BatchJob *job = itemToJob[id];
        batchJobManager->ShowLogFile(job);
    }
}

void JobsTree::CleanSuccessful()
{
    if (QMessageBox::question(this, "Confirm Delete Jobs", "Are you sure you want to delete all successful jobs?", QMessageBox::Yes | QMessageBox::No) != QMessageBox::Yes)
    {
        return;
    }
    batchJobManager->CleanSucc();
}

void JobsTree::CleanAll()
{
    if (QMessageBox::question(this, "Confirm Delete Jobs", "Are you sure you want to delete all jobs?", QMessageBox::Yes | QMessageBox::No) != QMessageBox::Yes)
    {
        return;
    }
    batchJobManager->CleanAll();
}

void JobsTree::OpenResults()
{
    QList<QTreeWidgetItem*> selected;
    GetSelections(selected);
    for (int i = 0; i < selected.size(); i++)
    {
        QTreeWidgetItem *id = selected[i];
        BatchJob *job = itemToJob[id];
        job->openResults();
    }
}

JobsView::JobsView(QWidget *parent)
    : QWidget(parent),
      listeningToBatchJobManager(false)
{

    jobsTree = new JobsTree(this);
    QVBoxLayout *vbox = new QVBoxLayout();
    vbox->setContentsMargins(0, 0, 0, 0);
    vbox->setSpacing(2);
    vbox->addWidget(jobsTree);
    setLayout(vbox);
}

JobsView::~JobsView()
{
}

void JobsView::jobAdded(BatchJob *job)
{
    addJobToTree(job);
    connect(job, SIGNAL(jobChanged(BatchJob*)),
            this, SLOT(jobChanged(BatchJob*)));
}

void JobsView::jobDeleted(BatchJob *job)
{
    QTreeWidgetItem *id = jobsTree->jobToItem[job];
    jobsTree->Delete(id);
    jobsTree->jobToItem.erase(job);
    jobsTree->itemToJob.erase(id);
}

void JobsView::jobChanged(BatchJob *job)
{
    QTreeWidgetItem *id = jobsTree->jobToItem[job];
    stylizeItem(id, job);
}

void JobsView::update(BatchJobManager *batchJobManager)
{
    if (!listeningToBatchJobManager)
    {
        listeningToBatchJobManager = true;
        connect(batchJobManager, SIGNAL(jobAdded(BatchJob*)),
                this, SLOT(jobAdded(BatchJob*)));
        connect(batchJobManager, SIGNAL(jobDeleted(BatchJob*)),
                this, SLOT(jobDeleted(BatchJob*)));
        jobsTree->batchJobManager = batchJobManager;
    }
    std::vector<BatchJob*> &jobList = *batchJobManager->GetJobList();
    std::vector<BatchJob*>::iterator jobIter = jobList.begin();
    while (jobIter != jobList.end())
    {
        BatchJob *job = *jobIter;
        ++jobIter;
        if (!jobInTree(job))
        {
            addJobToTree(job);
        }
    }
}

bool JobsView::jobInTree(BatchJob *job)
{
    return jobsTree->jobToItem.find(job) != jobsTree->jobToItem.end();
}

void JobsView::addJobToTree(BatchJob *job)
{
    QTreeWidgetItem *id = jobsTree->appendItem(jobsTree->rootId, "New Job", -1, -1);
    jobsTree->jobToItem[job] = id;
    jobsTree->itemToJob[id] = job;
    stylizeItem(id, job);
}

void JobsView::stylizeItem(QTreeWidgetItem *id, BatchJob *job)
{
    int image = 1;
    if (job->isRunning())
    {
        image = 2;
    }
    else if (job->isSuccess())
    {
        image = 3;
    }
    else
    {
        image = 4;
    }
    id->setIcon(0, jobsTree->GetIcon(image));
    QString jobText;
    if (!job->jobName().isEmpty())
    {
        jobText = QString("%1 (%2)").arg(job->jobName()).arg(job->jobId());
    }
    else
    {
        jobText = QString::number(job->jobId());
    }
    id->setText(0, jobText);
}


#include "JobsView.moc"
