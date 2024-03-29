#include "MainForm.h"

#include <QApplication>
#include <QClipboard>

#include "Render.h"

#include <boost/lexical_cast.hpp>
#include <numeric>

const int renderWidth = 920;
const int renderHeight = 920;
const int renderIterations = 999000;

MainForm::MainForm(QWidget* parent) :
	QWidget(parent)
{
	ui.setupUi(this);

	connect(ui.actionUpdate, &QAction::triggered, this, &MainForm::actionUpdateTriggered);
	connect(ui.actionRecalculate, &QAction::triggered, this, &MainForm::actionRecalculateTriggered);
	ui.splitter->setSizes({400, 920});

	connect(ui.actionUpdate, &QAction::triggered, this, &MainForm::actionUpdateTriggered);
	connect(ui.actionRecalculate, &QAction::triggered, this, &MainForm::actionRecalculateTriggered);

	connect(ui.valueV1button, &QPushButton::clicked, this, &MainForm::valueV1buttonClicked);
	connect(ui.valueV2button, &QPushButton::clicked, this, &MainForm::valueV2buttonClicked);
	connect(ui.valueV1button, &QPushButton::clicked, this, &MainForm::valueV1buttonClicked);
	connect(ui.valueV2button, &QPushButton::clicked, this, &MainForm::valueV2buttonClicked);

	connect(ui.checkBoxFindAllRoots, &QAbstractButton::toggled, this, &MainForm::checkBoxfindAllRootsToggled);

	ui.tableWidget->horizontalHeader()->setSectionResizeMode(QHeaderView::Stretch);
	ui.tableWidget->horizontalHeader()->setSectionsClickable(false);
	ui.tableWidget->verticalHeader()->setSectionResizeMode(QHeaderView::Fixed);
	connect(ui.tableWidget->verticalHeader(), &QHeaderView::sectionClicked, this, &MainForm::tableRowDoubleClicked);
	connect(ui.tableWidget, &QTableWidget::cellDoubleClicked, this, &MainForm::tableCellDoubleClicked);
	connect(ui.tableWidget, &QTableWidget::itemSelectionChanged, this, &MainForm::tableItemSelectionChanged);

	previewImage = QImage(renderWidth, renderHeight, QImage::Format_Grayscale8);

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

	bool warnB = false;
	if (nV1 != nV2 && nV1 != 2 && nV2 != 2) {
		// nB should be even for alternating patterns
		if (nB % 2 != 0) {
			// If nB is odd, we can try borrowing a *2 from nN
			if (nN %2 == 0 && nN > 2) {
				nA *= 2;
				nB *= 2;
				nN /= 2;
			}
			else {
				warnB = true;
			}
		}
	}
	ui.labelBvalue->setStyleSheet(warnB ? "background-color: red" : "");

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
	ui.tableWidget->clearSelection();
	ui.tableWidget->setRowCount(n);

	std::sort(roots.begin(), roots.end(), [](auto a, auto b){ return norm(a) >= norm(b); });

	for (int i = 0; i < n; i++) {
		if (roots[i].real() < 0)
			roots[i] = -roots[i];
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
		ui.tableWidget->selectRow(0);
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

void MainForm::tableItemSelectionChanged() {
	auto selection = ui.tableWidget->selectedItems();
	if (selection.empty())
		return;
	int row = selection.front()->row();
	if (row >= 0) {
		auto itemRe = ui.tableWidget->item(row, 0);
		auto itemIm = ui.tableWidget->item(row, 1);
		if (itemRe && itemIm) {
			complex_t root(boost::lexical_cast<double>(itemRe->text().toStdString()), boost::lexical_cast<double>(itemIm->text().toStdString()));
			Renderer renderer(renderWidth, renderHeight);
			Kleinian K(
				nA, nB, nN,
				ui.spinBoxV1->value(),
				ui.spinBoxV2->value()
			);
			Moebius<complex_t> camera{ 1, -1, 1, 1 }; // Halfplane to disc
			renderer.init(K, root, camera);
			renderer.clear();
			renderer.iterate(renderIterations);
			renderer.tonemap(previewImage);
			ui.labelPreview->setPixmap(QPixmap::fromImage(previewImage));
			return;
		}
	}
	clearPreview();
}

void MainForm::clearPreview() {
	QPixmap pixmap(renderWidth, renderHeight);
	pixmap.fill(Qt::black);
	ui.labelPreview->setPixmap(pixmap);
}

void MainForm::setPreview(QImage& image) {
	ui.labelPreview->setPixmap(QPixmap::fromImage(image));
}
