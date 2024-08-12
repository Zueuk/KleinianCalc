#pragma once

#include <QWidget>

#include "ui_MainForm.h"

#include "Moebius.hpp"

#include <complex>
using complex = std::complex<double>;

namespace Ui {
class MainForm;
}

class MainForm : public QWidget
{
	Q_OBJECT

	int nA, nB, nN;
	int nV1, nV2;
	QImage previewImage;

	Moebius<complex> Ma, Mb, Mc;

public:
	explicit MainForm(QWidget* parent = nullptr);

private:
	Ui::MainForm ui;

protected slots:
	void actionUpdateTriggered();
	void actionRecalculateTriggered();

	void valueUbuttonClicked();
	void valueV1buttonClicked();
	void valueV2buttonClicked();

	void copyXmlButtonClicked();

	void checkBoxfindAllRootsToggled(bool checked);

	void tableRowDoubleClicked(int index);
	void tableCellDoubleClicked(int row, int column);
	void tableItemSelectionChanged();

	void clearPreview();
	void setPreview(QImage& image);
};
