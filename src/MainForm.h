#pragma once

#include <QWidget>

#include "ui_MainForm.h"

namespace Ui {
class MainForm;
}

class MainForm : public QWidget
{
	Q_OBJECT

	int nA, nB, nN;

public:
	explicit MainForm(QWidget* parent = nullptr);

private:
	Ui::MainForm ui;

protected slots:
	void actionUpdateTriggered();
	void actionRecalculateTriggered();

	void valueV1buttonClicked();
	void valueV2buttonClicked();
	void checkBoxfindAllRootsToggled(bool checked);

	void tableRowDoubleClicked(int index);
	void tableCellDoubleClicked(int row, int column);
};
