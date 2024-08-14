#include "MainForm.h"

#include <QApplication>
#include <QClipboard>

#include "Render.h"

#include <boost/lexical_cast.hpp>

const int renderWidth = 920;
const int renderHeight = 920;
const int renderIterations = 999000;

MainForm::MainForm(QWidget* parent) :
	QWidget(parent)
{
	ui.setupUi(this);

	ui.splitter->setSizes({400, 920});
	ui.splitter->setStretchFactor(0, 1);

	connect(ui.actionUpdate, &QAction::triggered, this, &MainForm::actionUpdateTriggered);
	connect(ui.actionRecalculate, &QAction::triggered, this, &MainForm::actionRecalculateTriggered);

	ui.spinBoxU->hide();
	connect(ui.valueUbutton, &QPushButton::clicked, this, &MainForm::valueUbuttonClicked);
	connect(ui.valueV1button, &QPushButton::clicked, this, &MainForm::valueV1buttonClicked);
	connect(ui.valueV2button, &QPushButton::clicked, this, &MainForm::valueV2buttonClicked);

	connect(ui.checkBoxFindAllRoots, &QAbstractButton::toggled, this, &MainForm::checkBoxfindAllRootsToggled);

	connect(ui.buttonCopyXml, &QAbstractButton::clicked, this, &MainForm::copyXmlButtonClicked);

	ui.tableWidget->horizontalHeader()->setSectionResizeMode(QHeaderView::Stretch);
	ui.tableWidget->horizontalHeader()->setSectionsClickable(false);
	ui.tableWidget->verticalHeader()->setSectionResizeMode(QHeaderView::Fixed);
	connect(ui.tableWidget->verticalHeader(), &QHeaderView::sectionClicked, this, &MainForm::tableRowDoubleClicked);
	connect(ui.tableWidget, &QTableWidget::cellDoubleClicked, this, &MainForm::tableCellDoubleClicked);
	connect(ui.tableWidget, &QTableWidget::itemSelectionChanged, this, &MainForm::tableItemSelectionChanged);

	connect(ui.comboBoxViewMode, QOverload<int>::of(&QComboBox::currentIndexChanged), this, &MainForm::tableItemSelectionChanged);

	previewImage = QImage(renderWidth, renderHeight, QImage::Format_Grayscale8);

	Ma = { 0, 1, -1, 0 };
	Mb = { 0, 1, 1, 0 };
	Mc = { 0, 1, 1, 0 };

	actionUpdateTriggered();
}

void MainForm::actionUpdateTriggered() {
	K.a = ui.spinBoxA->value();
	K.b = ui.spinBoxB->value();
	K.n = ui.spinBoxN->value();

	// If GCD(a,b) > 1, we can move this multiplier into n
	int d = std::gcd(K.a, K.b);
	if (d > 1) {
		K.a /= d;
		K.b /= d;
		K.n *= d;
	}

	auto updateValueAndUi = [](int* value, QSpinBox* spinBox, QLabel* infLabel, bool infChecked) {
		if (infChecked) {
			*value = INT_MAX;
			spinBox->hide();
			infLabel->show();
		}
		else {
			*value = spinBox->value();
			infLabel->hide();
			spinBox->show();
		}
	};
	updateValueAndUi(&K.v1, ui.spinBoxV1, ui.labelInfV1, ui.checkBoxInfV1->isChecked());
	updateValueAndUi(&K.v2, ui.spinBoxV2, ui.labelInfV2, ui.checkBoxInfV2->isChecked());

	bool warnB = false;
	if (K.v1 != K.v2 && K.v1 > 2 && K.v2 > 2) {
		// b should be even for alternating patterns
		if (K.b % 2 != 0) {
			// If b is odd, we can try borrowing a *2 from n
			if (K.n > 2 && K.n % 2 == 0) {
				K.a *= 2;
				K.b *= 2;
				K.n /= 2;
			}
			else {
				warnB = true;
			}
		}
	}
	ui.labelBvalue->setStyleSheet(warnB ? "background-color: red" : "");

	ui.labelAvalue->setNum(K.a);
	ui.labelBvalue->setNum(K.b);
	ui.labelNvalue->setNum(K.n);

	// GCD(a,b)*n should be > 1
	ui.labelNvalue->setStyleSheet(K.n > 1 ? "" : "background-color: red");

	double v1 = offsetN<double>(K.v1);
	double v2 = offsetN<double>(K.v2);
	Mb.d.imag(v1);
	Mc.d.imag(-v2);
	QString v1string = QString::number(v1, 'g', 8);
	QString v2string = QString::number(v2, 'g', 8);
	ui.valueV1button->setText(v1string);
	ui.valueV2button->setText(v2string);

	ui.label_b22->setText("0\n" + v1string);
	ui.label_B11->setText("0\n" + v1string);

	ui.label_c22->setText("0\n-" + v2string);
	ui.label_C11->setText("0\n-" + v2string);

	if (K.a + K.b < 42) {
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

	int method = ui.radioBtnMethod0->isChecked() ? 0 : 1;
	bool findAllRoots = ui.checkBoxFindAllRoots->isChecked();
	auto roots = K.solve(method, findAllRoots);

	int n = (int)roots.size();
	ui.tableWidget->clearSelection();
	ui.tableWidget->setRowCount(n);

	std::sort(roots.begin(), roots.end(), [](auto a, auto b){ return norm(a) >= norm(b); });

	for (int i = 0; i < n; ++i) {
		// Real part must be positive for the pattern to be on the right side
		if (roots[i].real() < 0)
			roots[i] = -roots[i];
		// Very small imaginary value must be just approximation error
		if (fabs(roots[i].imag()) < 1e-10)
			roots[i].imag(0);

		auto root = complex(roots[i]);
		auto itemRe = new QTableWidgetItem(QString::fromStdString(boost::lexical_cast<std::string>(root.real())));
		auto itemIm = new QTableWidgetItem(QString::fromStdString(boost::lexical_cast<std::string>(root.imag())));
		itemRe->setData(Qt::UserRole, root.real());
		itemIm->setData(Qt::UserRole, root.imag());
		if (abs(root) < 1 || abs(root) > 2) {
			itemRe->setForeground(Qt::darkGray);
			itemIm->setForeground(Qt::darkGray);
		}
		ui.tableWidget->setItem(i, 0, itemRe);
		ui.tableWidget->setItem(i, 1, itemIm);
	}

	if (n != 0) {
		ui.tableWidget->selectRow(0);
		ui.buttonCopyXml->setEnabled(true);
	}
	else {
		clearPreview();
		ui.buttonCopyXml->setEnabled(false);
	}

	QApplication::restoreOverrideCursor();
}

static QString toParamString(const QString& name, const complex& z) {
	return "Re_" + name + "=\"" + QString::number(z.real(), 'g', 16) + "\" "
		"Im_" + name + "=\"" + QString::number(z.imag(), 'g', 16) + "\" ";
}

static QString toXformString(const Moebius<complex>& M) {
	return "<xform weight=\"1\" color=\"0\" mobius=\"1\" coefs=\"1 0 0 1 0 0\" " +
		toParamString("A", M.a) +
		toParamString("B", M.b) +
		toParamString("C", M.c) +
		toParamString("D", M.d) +
		"/>\n";
}

void MainForm::copyXmlButtonClicked() {
	QString viewScale, viewXform;
	if (ui.comboBoxViewMode->currentIndex() == 0) {
		viewScale = "512";
		viewXform = "<finalxform color=\"0\" symmetry=\"1\" mobius=\"1\" coefs=\"1 0 0 1 0 0\""
					" Re_A=\"1\" Im_A=\"0\" Re_B=\"-1\" Im_B=\"0\" Re_C=\"1\" Im_C=\"0\" Re_D=\"1\" Im_D=\"0\"/>\n";
	} else {
		viewScale = "256";
		viewXform = "";
	}

	QApplication::clipboard()->setText(
		"<flame name=\"Kleinian\""
		" size=\"1024 1024\" center=\"0 0\" scale=\"" + viewScale + "\""
		" oversample=\"1\" filter=\"0.5\" quality=\"50\""
		" background=\"0 0 0\" brightness=\"4\" gamma=\"4\" gamma_threshold=\"0.04\">\n" +
			toXformString(Ma) +
			toXformString(Ma.inverse()) +
			toXformString(Mb) +
			toXformString(Mb.inverse()) +
			toXformString(Mc) +
			toXformString(Mc.inverse()) +
			viewXform +
			"<palette count=\"1\" format=\"RGB\">"
				"FFFFFF"
			"</palette>\n"
		"</flame>"
	);
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

void MainForm::valueUbuttonClicked() {
	double u = 2.0;
	QApplication::clipboard()->setText(QString::fromStdString(boost::lexical_cast<std::string>(u)));
}

void MainForm::valueV1buttonClicked() {
	double v1 = offsetN<double>(K.v1);
	QApplication::clipboard()->setText(QString::fromStdString(boost::lexical_cast<std::string>(v1)));
}

void MainForm::valueV2buttonClicked() {
	double v2 = offsetN<double>(K.v2);
	QApplication::clipboard()->setText(QString::fromStdString(boost::lexical_cast<std::string>(v2)));
}

void MainForm::tableItemSelectionChanged() {
	auto selection = ui.tableWidget->selectedItems();
	if (selection.isEmpty()) {
		ui.label_a22->setText("...");
		ui.label_A11->setText("...");
		return;
	}

	int row = selection.front()->row();
	if (row >= 0) {
		auto itemRe = ui.tableWidget->item(row, 0);
		auto itemIm = ui.tableWidget->item(row, 1);
		if (itemRe && itemIm) {
			double re = itemRe->data(Qt::UserRole).value<double>();
			double im = itemIm->data(Qt::UserRole).value<double>();
			Ma.d = complex(re, im);
			QString rootString = QString::number(re, 'g', 8) + "\n" + QString::number(im, 'g', 8);
			ui.label_a22->setText(rootString);
			ui.label_A11->setText(rootString);

			Renderer::complex root(re, im);
			Renderer renderer(renderWidth, renderHeight);
			static const Moebius<Renderer::complex> views[] = {
				{ 1, -1, 1, 1 }, // Halfplane to disc
				{ 0.5, 0, 0, 1 }, // Zoom out
			};
			int viewIndex = qBound(0, ui.comboBoxViewMode->currentIndex(), 1);
			renderer.init(K, root, views[viewIndex]);
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
