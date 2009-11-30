#ifndef MIQTREEWIDGET_H
#define MIQTREEWIDGET_H

#include <QTreeWidget>
#include <QList>
#include <QIcon>

#include <vector>
#include <string>

class TreeData;

class MIQTreeWidget : public QTreeWidget
{
    Q_OBJECT
public:
    MIQTreeWidget(QWidget *parent) : QTreeWidget(parent)
    {
    }
    ~MIQTreeWidget();

    void GetSelections(QList<QTreeWidgetItem*> &selected);
    TreeData *GetItemData(QTreeWidgetItem *item);
    void AssignImageList(std::vector<QIcon> &imageList);

    void Delete(QTreeWidgetItem *item);
    void DeleteChildren(QTreeWidgetItem *item);

    QTreeWidgetItem *appendItem(QTreeWidgetItem *parent,
                                const std::string &name,
                                int image = -1,
                                int selImage = -1,
                                TreeData *data = 0);

    QTreeWidgetItem *insertItem(QTreeWidgetItem *parent,
                                QTreeWidgetItem *previous,
                                const std::string &name,
                                int image = -1,
                                int selImage = -1,
                                TreeData *data = 0);

    QTreeWidgetItem *prependItem(QTreeWidgetItem *parent,
                                 const std::string &name,
                                 int image = -1,
                                 int selImage = -1,
                                 TreeData *data = 0);

    QIcon&GetIcon(unsigned int image_num);

protected:

    std::vector<QIcon> _imageList;

private:
    void itemInit(QTreeWidgetItem *item, const std::string &name, int image, TreeData *data);
};


// this is required so that TreeData can be stored in a QVariant
Q_DECLARE_METATYPE(TreeData*);

#endif // ifndef MIQTREEWIDGET_H
