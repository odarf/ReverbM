#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include "processing.h"
#include "analysis.h"
#include <QMainWindow>
#include <QtWidgets>
#include <QMediaPlayer>

QT_BEGIN_NAMESPACE
namespace Ui { class MainWindow; }
QT_END_NAMESPACE

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    MainWindow(QWidget *parent = nullptr);
    float delay, decay, mixPercent;
    ~MainWindow();


private slots:
    void on_pushButton_clicked();

    void on_delaySlider_valueChanged(int value);

    void on_decaySlider_valueChanged(int value);

    void on_dwSlider_valueChanged(int value);


private:
    Ui::MainWindow *ui;
};
#endif // MAINWINDOW_H
