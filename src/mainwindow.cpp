#include "mainwindow.h"
#include "ui_mainwindow.h"
#include "modeling.h"
#include "analysis.h"
#include "inou.h"
#include "qmath.h"

#include <iostream>
#include <QLogValueAxis>
#include <QLineSeries>
#include <QValueAxis>
#include <QChart>
#include <QChartView>
#include <QtCharts>
#include <QtWidgets>
#include <QtWidgets/QHBoxLayout>
#include <QHBoxLayout>
#include <cstdlib>

#include <math.h>
#include <time.h>
#include <random>

#include <fstream>

using namespace std;

static const char* INPUT_WAV = "../soundInput.wav";
static const char* OUTPUT_WAV = "../Reverbed.wav";

MainWindow::MainWindow(QWidget *parent): QMainWindow(parent), ui(new Ui::MainWindow) {
    ui->setupUi(this);
}

MainWindow::~MainWindow() { delete ui; }

void MainWindow::on_pushButton_clicked(){
    modeling  modeling;
    processing processing;
    analysis analysis;
    QVector<double> x, reverb, wavReverbed;
    double max, min;

    QVector<double> wav = inou().loadWave(INPUT_WAV);

    for (int i(0); i<wav.length(); i++) { x.append(i); }

    ui->plotInputWave->addGraph();
    ui->plotInputWave->graph(0)->setData(x, wav);
    ui->plotInputWave->xAxis->setRange(0, wav.length());
    max = analysis.maxValue(wav);
    min = analysis.minValue(wav);
    ui->plotInputWave->yAxis->setRange(min, max);
    ui->plotInputWave->replot();

    int delay = 500;
    int delaySamples = (int)((float)delay * 22.050f);
    float decay = 0.5f;

    for (int i(0); i<wav.length(); i++) {
        wav[i+delaySamples] += (short)((float)wav[i] * decay);
    }

    ui->plotReverbedWave->addGraph();
    ui->plotReverbedWave->graph(0)->setData(x, wav);
    ui->plotReverbedWave->xAxis->setRange(0, wav.length());
    max = analysis.maxValue(wav);
    min = analysis.minValue(wav);
    ui->plotReverbedWave->yAxis->setRange(min, max);
    ui->plotReverbedWave->replot();

    inou().exportWave(wav, wav.length(), OUTPUT_WAV, 1);
}

