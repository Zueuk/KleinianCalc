#pragma once

#include <QWidget>

#include "ui_MainForm.h"

#include "Kleinian.hpp"
#include "Renderer.h"

#include <complex>
using complex = std::complex<double>;

namespace Ui {
class MainForm;
}

class MainForm : public QWidget
{
	Q_OBJECT

	Kleinian K;

	Renderer renderer;
	QImage previewImage;

	Moebius<complex> Ma, Mb, Mc, Mseq1, Mseq2;
	double offsetH, offsetV1, offsetV2;

	Moebius<complex> sequenceMatrix(const std::string& sequence);

public:
	explicit MainForm(QWidget* parent = nullptr);

private:
	Ui::MainForm ui;

protected slots:
	void actionUpdateTriggered();
	void actionRecalculateTriggered();

	void valueHbuttonClicked();
	void valueV1buttonClicked();
	void valueV2buttonClicked();

	void copyImageButtonClicked();
	void copyXmlButtonClicked();

	void checkBoxfindAllRootsToggled(bool checked);

	void tableRowDoubleClicked(int index);
	void tableCellDoubleClicked(int row, int column);
	void tableItemSelectionChanged();

	void updateMseq1(const QString& str);
	void updateMseq2(const QString& str);
	void lineEditSeq1Edited(const QString& str);
	void lineEditSeq2Edited(const QString& str);

	void clearPreview();
	void renderPreview();
};
