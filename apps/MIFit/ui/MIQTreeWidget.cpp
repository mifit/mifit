#include "MIQTreeWidget.h"
#include "TreeData.h"

using namespace boost::signals;

MIQTreeWidget::~MIQTreeWidget() {
  SignalConnectionList::iterator iter;
  for (iter = signalConnections.begin(); iter != signalConnections.end(); ++iter) {
    iter->disconnect();
  }
}


void MIQTreeWidget::AssignImageList(std::vector<QIcon> &imageList) {
  _imageList=imageList;
  for (size_t i = 0; i < _imageList.size(); ++i) {

    // override automatically-generated icon modes (selected, etc) with the default one
    QPixmap p=_imageList[i].pixmap(16,16);
    _imageList[i].addPixmap(p, QIcon::Selected, QIcon::On);
    _imageList[i].addPixmap(p, QIcon::Selected, QIcon::Off);
  }
}


void MIQTreeWidget::GetSelections(QList<QTreeWidgetItem *> &selected) {
  selected=selectedItems();
}

TreeData* MIQTreeWidget::GetItemData(QTreeWidgetItem *item) {
  if (!item)
    return 0;
  return item->data(0,Qt::UserRole).value<TreeData*>();
}

void MIQTreeWidget::Delete(QTreeWidgetItem *item) {
  //FIXME is this correct?
  delete item;
}

void MIQTreeWidget::DeleteChildren(QTreeWidgetItem *item) {
  QList<QTreeWidgetItem*> items=item->takeChildren();
  for (int i=0; i < items.size(); ++i) {
    Delete(items[i]);
  }
}

void MIQTreeWidget::itemInit(QTreeWidgetItem *item, const std::string &name, int image, TreeData *data)
{
  if (image != -1)
    item->setIcon(0,_imageList[image]);
  item->setText(0,name.c_str());

  QVariant v;
  v.setValue<TreeData*>(data);
  item->setData(0,Qt::UserRole,v);
  if (data)
    data->SetId(item);
}

QTreeWidgetItem *MIQTreeWidget::appendItem(QTreeWidgetItem *parent,
                                           const std::string &name,
                                           int image,
                                           int,
                                           TreeData *data) {

  QTreeWidgetItem *item = new QTreeWidgetItem(parent);
  itemInit(item,name,image,data);
  return item;
}

QTreeWidgetItem *MIQTreeWidget::insertItem(QTreeWidgetItem *parent,
                                           QTreeWidgetItem *previous,
                                           const std::string &name,
                                           int image,
                                           int ,
                                           TreeData *data) {
  QTreeWidgetItem *item = new QTreeWidgetItem(parent,previous);
  itemInit(item,name,image,data);
  return item;
}

QTreeWidgetItem *MIQTreeWidget::prependItem(QTreeWidgetItem *parent,
                                            const std::string &name,
                                            int image,
                                            int,
                                            TreeData *data) {
  QTreeWidgetItem *item = new QTreeWidgetItem(parent);
  itemInit(item,name,image,data);
  return item;
}

void MIQTreeWidget::addSignalConnection(connection c) {
  c.set_controlling(false);
  signalConnections.push_back(c);
}

QIcon &MIQTreeWidget::GetIcon(unsigned int id) {
  static QIcon dummy;
  if (id > _imageList.size())
    return dummy;
  return _imageList[id];
}

