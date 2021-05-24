#include "MainForm.h"

#include <QApplication>
#include <QClipboard>

#include "Kleinian.h"

#include <boost/lexical_cast.hpp>
#include <numeric>

MainForm::MainForm(QWidget* parent) :
	QWidget(parent)
{
	ui.setupUi(this);

	connect(ui.actionUpdate, &QAction::triggered, this, &MainForm::actionUpdateTriggered);
	connect(ui.actionRecalculate, &QAction::triggered, this, &MainForm::actionRecalculateTriggered);

	connect(ui.valueV1button, &QPushButton::clicked, this, &MainForm::valueV1buttonClicked);
	connect(ui.valueV2button, &QPushButton::clicked, this, &MainForm::valueV2buttonClicked);

	connect(ui.checkBoxFindAllRoots, &QAbstractButton::toggled, this, &MainForm::checkBoxfindAllRootsToggled);

	ui.tableWidget->horizontalHeader()->setSectionResizeMode(QHeaderView::Stretch);
	ui.tableWidget->horizontalHeader()->setSectionsClickable(false);
	ui.tableWidget->verticalHeader()->setSectionResizeMode(QHeaderView::Fixed);
	connect(ui.tableWidget->verticalHeader(), &QHeaderView::sectionClicked, this, &MainForm::tableRowDoubleClicked);
	connect(ui.tableWidget, &QTableWidget::cellDoubleClicked, this, &MainForm::tableCellDoubleClicked);

	actionUpdateTriggered();
}

void MainForm::actionUpdateTriggered() {
	nA = ui.spinBoxA->value();
	nB = ui.spinBoxB->value();
	nN = ui.spinBoxN->value();

	// If GCD(a,b) > 1, we can move this multiplier into n
	int d = std::gcd(nA, nB);
	if (d > 1) {
		nA /= d;
		nB /= d;
		nN *= d;
	}

	int nV1 = ui.spinBoxV1->value();
	int nV2 = ui.spinBoxV2->value();

	ui.labelAvalue->setNum(nA);
	ui.labelBvalue->setNum(nB);
	ui.labelNvalue->setNum(nN);

	// GCD(a,b)*n should be > 1
	ui.labelNvalue->setStyleSheet(nN > 1 ? "" : "background-color: red");

	double v1 = offsetN<double>(nV1);
	double v2 = offsetN<double>(nV2);
	QString v1string = QString::number(v1, 'g', 8);
	QString v2string = QString::number(v2, 'g', 8);
	ui.valueV1button->setText(v1string);
	ui.valueV2button->setText(v2string);

	ui.label_b22->setText("0\n" + v1string);
	ui.label_B11->setText("0\n" + v1string);

	ui.label_c22->setText("0\n-" + v2string);
	ui.label_C11->setText("0\n-" + v2string);

	if (nA + nB < 42) {
		// If the polynomial is not too big, solve it right away
		ui.buttonRecalc->hide();
		actionRecalculateTriggered();
	}
	else {
		// Otherwise, let the user press "Solve" when ready
		ui.tableWidget->setRowCount(0);
		ui.buttonRecalc->show();
	}
}

void MainForm::actionRecalculateTriggered() {
	// Solving can take some time
	QApplication::setOverrideCursor(QCursor(Qt::WaitCursor));

	Kleinian K(
		nA, nB, nN,
		ui.spinBoxV1->value(),
		ui.spinBoxV2->value()
	);

	int method = ui.radioBtnMethod0->isChecked() ? 0 : 1;
	bool findAllRoots = ui.checkBoxFindAllRoots->isChecked();
	auto roots = K.solve(method, findAllRoots);

	int n = (int)roots.size();
	ui.tableWidget->setRowCount(n);

	std::sort(roots.begin(), roots.end(), [](auto a, auto b){ return abs(a) >= abs(b); });

	for (int i = 0; i < n; i++) {
		auto root = std::complex<double>(roots[i]);
		auto itemRe = new QTableWidgetItem(QString::fromStdString(boost::lexical_cast<std::string>(root.real())));
		auto itemIm = new QTableWidgetItem(QString::fromStdString(boost::lexical_cast<std::string>(root.imag())));
		if (abs(root) < 1 || abs(root) > 2) {
			itemRe->setForeground(Qt::darkGray);
			itemIm->setForeground(Qt::darkGray);
		}
		ui.tableWidget->setItem(i, 0, itemRe);
		ui.tableWidget->setItem(i, 1, itemIm);
	}

	if (roots.size() != 0) {
		auto root = std::complex<double>(roots.front());
		QString rootString = QString::number(root.real(), 'g', 8) + "\n" + QString::number(root.imag(), 'g', 8);
		ui.label_a22->setText(rootString);
		ui.label_A11->setText(rootString);
	}
	else {
		ui.label_a22->setText("...");
		ui.label_A11->setText("...");
	}

	QApplication::restoreOverrideCursor();
}

void MainForm::checkBoxfindAllRootsToggled(bool checked) {
	if (checked && ui.tableWidget->rowCount() < 2) {
		actionUpdateTriggered();
	}
}

void MainForm::tableRowDoubleClicked(int row) {
	auto item = ui.tableWidget->item(row, 0);
	if (item) {
		QString text = item->text() + " " + ui.tableWidget->item(row, 1)->text();
		QApplication::clipboard()->setText(text);
	}
}

void MainForm::tableCellDoubleClicked(int row, int column) {
	auto item = ui.tableWidget->item(row, column);
	if (item) {
		QString text = item->text();
		QApplication::clipboard()->setText(text);
	}
}

void MainForm::valueV1buttonClicked() {
	double v1 = offsetN<double>(ui.spinBoxV1->value());
	QApplication::clipboard()->setText(QString::fromStdString(boost::lexical_cast<std::string>(v1)));
}

void MainForm::valueV2buttonClicked() {
	double v2 = offsetN<double>(ui.spinBoxV2->value());
	QApplication::clipboard()->setText(QString::fromStdString(boost::lexical_cast<std::string>(v2)));
}
