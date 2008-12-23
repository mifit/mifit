#include "EnvironmentPreferences.h"
#include "corelib.h"
#include "uilib.h"

#include <QFileDialog>

EnvironmentPreferences::EnvironmentPreferences(QWidget *parent) : QWidget(parent) {
  setupUi(this);
  Application* app = Application::instance();

  crystalDataLineEdit->setText(app->CrystalData.c_str());
  customDictionaryLineEdit->setText(app->XFitDictSetting.c_str());
  shelxHomeLineEdit->setText(app->ShelxHome.c_str());
  htmlBrowserLineEdit->setText(app->HTMLBrowser.c_str());
  smilesDatabaseLineEdit->setText(app->SmilesDbCommand.c_str());
  checkpointDirectoryLineEdit->setText(app->checkpointDirectory.c_str());
}

void EnvironmentPreferences::savePreferences() {
  Application* app = Application::instance();
  app->CrystalData = crystalDataLineEdit->text().toStdString();
  app->XFitDictSetting = customDictionaryLineEdit->text().toStdString();
  app->ShelxHome = shelxHomeLineEdit->text().toStdString();
  app->HTMLBrowser = htmlBrowserLineEdit->text().toStdString();
  app->SmilesDbCommand = smilesDatabaseLineEdit->text().toStdString();
  app->checkpointDirectory = checkpointDirectoryLineEdit->text().toStdString();
}

void EnvironmentPreferences::on_crystalDataButton_pressed() {
  crystalDataLineEdit->setText(QFileDialog::getExistingDirectory(this, "", crystalDataLineEdit->text()));
}
void EnvironmentPreferences::on_customDictionaryButton_pressed() {
  customDictionaryLineEdit->setText(QFileDialog::getOpenFileName(this, "", customDictionaryLineEdit->text()));
}
void EnvironmentPreferences::on_shelxHomeButton_pressed() {
  shelxHomeLineEdit->setText(QFileDialog::getExistingDirectory(this, "", shelxHomeLineEdit->text()));
}
void EnvironmentPreferences::on_htmlBrowserButton_pressed() {
  htmlBrowserLineEdit->setText(QFileDialog::getOpenFileName(this, "", htmlBrowserLineEdit->text()));
}
void EnvironmentPreferences::on_checkpointDirectoryButton_pressed() {
  checkpointDirectoryLineEdit->setText(QFileDialog::getExistingDirectory(this, "", checkpointDirectoryLineEdit->text()));
}
