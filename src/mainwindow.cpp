#include "mainwindow.h"
#include "ui_mainwindow.h"
#include "analysis.h"
#include "reverberator.h"
#include "inou.h"

using namespace std;

static const char* INPUT_WAV = "../soundInput.wav";
static const char* OUTPUT_WAV = "../Reverbed.wav";

template<typename T> vector<T> arrange(T start, T stop, T step = 1) {
    vector<T> values;
    for (T value = start; value < stop; value += step)
        values.push_back(value);
    return values;
}

MainWindow::MainWindow(QWidget *parent): QMainWindow(parent), ui(new Ui::MainWindow) {
    ui->setupUi(this);
    ui->sbDecay->setValue(ui->decaySlider->value()/100.0f); //hack
    connect(ui->sbDelay, SIGNAL(valueChanged(int)), ui->delaySlider, SLOT(setValue(int)));
    connect(ui->sbDW, SIGNAL(valueChanged(int)), ui->dwSlider, SLOT(setValue(int)));

}

MainWindow::~MainWindow() { delete ui; }

void MainWindow::on_pushButton_clicked(){
    analysis analysis;
    QVector<double> x;
    double max, min;
    double dt = 1/22.050;

    QVector<double> wav = inou().loadWave(INPUT_WAV);

    for (int i(0); i<wav.length(); i++) { x.append(i); }

    ui->plotInputWave->addGraph();
    ui->plotInputWave->graph(0)->setData(x, wav);
    ui->plotInputWave->xAxis->setRange(0, wav.length());
    max = analysis.maxValue(wav);
    min = analysis.minValue(wav);
    ui->plotInputWave->yAxis->setRange(min, max);
    ui->plotInputWave->replot();

    QVector<double> reverbed = Reverberator().reverb(wav, this->delay, this->decay, this->mixPercent);
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
    inou().exportWave(reverbed, reverbed.length(), OUTPUT_WAV, 1);
    inou().exportWave(wav, wav.length(), "../test.wav", 1);

    wav = analysis.fourierAmplitude(wav);
    vector<double> temp = arrange<double>(0, (double)wav.length() * dt, dt);
    QVector<double> arranged;

    for (int i = 0; i<temp.size(); i++) {
        arranged.append(temp[i]);
    }

    ui->plotInputFour->addGraph();
    ui->plotInputFour->graph(0)->setData(arranged, wav);
    ui->plotInputFour->xAxis->setRange(0, wav.length());
    ui->plotInputFour->setInteraction(QCP::iRangeZoom,true);
    ui->plotInputFour->setInteraction(QCP::iRangeDrag,true);
    max = analysis.maxValue(wav);
    min = analysis.minValue(wav);
    ui->plotInputFour->yAxis->setRange(min, max);
    ui->plotInputFour->replot();

    reverbed = analysis.fourierAmplitude(reverbed);
    temp = arrange<double>(0, (double)reverbed.length() * dt, dt);
    arranged.clear();

    for (int i = 0; i<temp.size(); i++) {
        arranged.append(temp[i]);
    }

    ui->plotReverbedFour->addGraph();
    ui->plotReverbedFour->graph(0)->setData(arranged, reverbed);
    ui->plotReverbedFour->xAxis->setRange(0, reverbed.length());
    ui->plotReverbedFour->setInteraction(QCP::iRangeZoom,true);
    ui->plotReverbedFour->setInteraction(QCP::iRangeDrag,true);
    max = analysis.maxValue(reverbed);
    min = analysis.minValue(reverbed);
    ui->plotReverbedFour->yAxis->setRange(min, max);
    ui->plotReverbedFour->replot();

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

