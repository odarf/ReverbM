#include "mainwindow.h"
#include "ui_mainwindow.h"
#include "modeling.h"
#include "analysis.h"
#include "reverberator.h"
#include "inou.h"
#include "qmath.h"

#include <iostream>
#include <QLogValueAxis>
#include <QLineSeries>
#include <QValueAxis>
#include <QChart>

#include <QSound>
#include <QMediaPlayer>
#include <QMediaPlaylist>

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
    ui->sbDecay->setValue(ui->decaySlider->value()/100.0f); //hack
    connect(ui->sbDelay, SIGNAL(valueChanged(int)), ui->delaySlider, SLOT(setValue(int)));
    connect(ui->sbDW, SIGNAL(valueChanged(int)), ui->dwSlider, SLOT(setValue(int)));

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

    QVector<double> reverbed = Reverberator().reverb(wav, this->delay, this->decay, this->mixPercent);
    //QVector<double> reverbed = Reverberator().reverb(wav, 78.9f, 0.45f, 50);
    double maxelem = analysis.maxValue(wav);
    for (int i = 0; i<reverbed.length(); i++) {
        reverbed[i] *= maxelem;
    }

    ui->plotReverbedWave->addGraph();
    ui->plotReverbedWave->graph(0)->setData(x, reverbed);
    ui->plotReverbedWave->xAxis->setRange(0, reverbed.length());
    max = analysis.maxValue(reverbed);
    min = analysis.minValue(reverbed);
    ui->plotReverbedWave->yAxis->setRange(min, max);
    ui->plotReverbedWave->replot();
    inou().exportWave(reverbed, reverbed.length(), "../combRev.wav", 1);
    clog << reverbed.length() << endl;
    clog << wav.length() << endl;
    inou().exportWave(wav, wav.length(), OUTPUT_WAV, 1);
}


void MainWindow::on_delaySlider_valueChanged(int value){
    delay = value/1.0f;
    ui->sbDelay->setValue(value);
}


void MainWindow::on_decaySlider_valueChanged(int value){
    decay = value/100.0f;
    ui->sbDecay->setValue(decay);
}


void MainWindow::on_dwSlider_valueChanged(int value){
    mixPercent = value;
    ui->sbDW->setValue(value);
}

void MainWindow::on_sbDecay_valueChanged(double arg1)
{
    ui->decaySlider->setValue((int)(arg1 * 100));
}

