#include "MainForm.h"

#include <QApplication>
#include <QClipboard>
#include <QRegExpValidator>

#include <boost/lexical_cast.hpp>

const int renderWidth = 920;
const int renderHeight = 920;
const int renderIterations = 2000000;

MainForm::MainForm(QWidget* parent)
	: QWidget(parent)
	, renderer(renderWidth, renderHeight)
	, previewImage(renderWidth, renderHeight, QImage::Format_Grayscale8)
{
	ui.setupUi(this);

	ui.splitter->setSizes({400, 920});
	ui.splitter->setStretchFactor(0, 1);

	connect(ui.actionUpdate, &QAction::triggered, this, &MainForm::actionUpdateTriggered);
	connect(ui.actionRecalculate, &QAction::triggered, this, &MainForm::actionRecalculateTriggered);

	connect(ui.valueHbutton, &QPushButton::clicked, this, &MainForm::valueHbuttonClicked);
	connect(ui.valueV1button, &QPushButton::clicked, this, &MainForm::valueV1buttonClicked);
	connect(ui.valueV2button, &QPushButton::clicked, this, &MainForm::valueV2buttonClicked);

	connect(ui.checkBoxFindAllRoots, &QAbstractButton::toggled, this, &MainForm::checkBoxfindAllRootsToggled);

	connect(ui.buttonCopyImage, &QAbstractButton::clicked, this, &MainForm::copyImageButtonClicked);
	connect(ui.buttonCopyXml, &QAbstractButton::clicked, this, &MainForm::copyXmlButtonClicked);

	ui.tableWidget->horizontalHeader()->setSectionResizeMode(QHeaderView::Stretch);
	ui.tableWidget->horizontalHeader()->setSectionsClickable(false);
	ui.tableWidget->verticalHeader()->setSectionResizeMode(QHeaderView::Fixed);
	connect(ui.tableWidget->verticalHeader(), &QHeaderView::sectionClicked, this, &MainForm::tableRowDoubleClicked);
	connect(ui.tableWidget, &QTableWidget::cellDoubleClicked, this, &MainForm::tableCellDoubleClicked);
	connect(ui.tableWidget, &QTableWidget::itemSelectionChanged, this, &MainForm::tableItemSelectionChanged);

	connect(ui.comboBoxCalcMode, QOverload<int>::of(&QComboBox::currentIndexChanged), this, &MainForm::tableItemSelectionChanged);
	connect(ui.comboBoxViewMode, QOverload<int>::of(&QComboBox::currentIndexChanged), this, &MainForm::tableItemSelectionChanged);

	connect(ui.groupBox_a, &QGroupBox::toggled, this, &MainForm::renderPreview);
	connect(ui.groupBox_A, &QGroupBox::toggled, this, &MainForm::renderPreview);
	connect(ui.groupBox_b, &QGroupBox::toggled, this, &MainForm::renderPreview);
	connect(ui.groupBox_B, &QGroupBox::toggled, this, &MainForm::renderPreview);
	connect(ui.groupBox_c, &QGroupBox::toggled, this, &MainForm::renderPreview);
	connect(ui.groupBox_C, &QGroupBox::toggled, this, &MainForm::renderPreview);
	connect(ui.groupBox_x, &QGroupBox::toggled, this, &MainForm::renderPreview);
	connect(ui.groupBox_X, &QGroupBox::toggled, this, &MainForm::renderPreview);

	auto abcValidator = new QRegExpValidator(QRegExp("[abcABC]*"), this);
	ui.lineEditSeq1->setValidator(abcValidator);
	ui.lineEditSeq2->setValidator(abcValidator);
	connect(ui.lineEditSeq1, &QLineEdit::textEdited, this, &MainForm::lineEditSeq1Edited);
	connect(ui.lineEditSeq2, &QLineEdit::textEdited, this, &MainForm::lineEditSeq2Edited);

	actionUpdateTriggered();
}

Moebius<complex> MainForm::sequenceMatrix(const std::string& sequence) {
	Moebius<complex> M;
	for (char ch : sequence) {
		switch (ch) {
			case 'a':
				M = M * Ma;
				break;
			case 'b':
				M = M * Mb;
				break;
			case 'c':
				M = M * Mc;
				break;
			case 'A':
				M = M * Ma.inverse();
				break;
			case 'B':
				M = M * Mb.inverse();
				break;
			case 'C':
				M = M * Mc.inverse();
				break;
		}
	}
	if (abs(M.c) >= FLT_EPSILON)
		M.divideBy(M.c);
	else
		M.divideBy(M.a);
	cleanMatrix(M);
	return M;
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

	auto updateValueAndUi = [](int* value, QSpinBox* spinBox, QLabel* infLabel, const QCheckBox* checkBox) {
		if (checkBox->isChecked()) {
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
	updateValueAndUi(&K.h, ui.spinBoxH, ui.labelInfH, ui.checkBoxInfH);
	updateValueAndUi(&K.v1, ui.spinBoxV1, ui.labelInfV1, ui.checkBoxInfV1);
	updateValueAndUi(&K.v2, ui.spinBoxV2, ui.labelInfV2, ui.checkBoxInfV2);

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

	if (isInf(K.h)) {
		offsetH = offsetN<double>(K.h);
		offsetV1 = offsetN<double>(K.v1);
		offsetV2 = offsetN<double>(K.v2);
	}
	else {
		TilngParameters<double> tiling(K.h, K.v1, K.v2);
		offsetH= tiling.ph;
		offsetV1 = tiling.qh;
		offsetV2 = tiling.rh;
	}
	ui.valueHbutton->setText(QString::number(offsetH, 'g', 8));
	ui.valueV1button->setText(QString::number(offsetV1, 'g', 8));
	ui.valueV2button->setText(QString::number(offsetV2, 'g', 8));

	int na = (K.b > 0) ? K.a / K.b : K.a;
	QString seq = QString(na + 1, 'a');
	ui.lineEditSeq1->setText(seq);
	ui.lineEditSeq2->setText(seq.toUpper());

	if (K.a + K.b < 42) {
		// If the polynomial is not too big, solve it right away
		ui.buttonRecalc->hide();
		actionRecalculateTriggered();
	}
	else {
		// Otherwise, let the user click "Solve" when ready
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

	ui.tableWidget->clearContents();
	ui.tableWidget->setRowCount(n);

	for (int i = 0; i < n; ++i) {
		auto root = complex(roots[i]);
		if (isInf(K.h)) {
			// Real part must be positive for the pattern to be on the right side
			if (root.real() < 0)
				root.real(-root.real());
			// Very small imaginary value must be just approximation error
			if (abs(root.imag()) < FLT_EPSILON)
				root.imag(0);
		}
		auto itemRe = new QTableWidgetItem(QString::fromStdString(boost::lexical_cast<std::string>(root.real())));
		auto itemIm = new QTableWidgetItem(QString::fromStdString(boost::lexical_cast<std::string>(root.imag())));
		itemRe->setData(Qt::UserRole, root.real());
		itemIm->setData(Qt::UserRole, root.imag());
		ui.tableWidget->setItem(i, 0, itemRe);
		ui.tableWidget->setItem(i, 1, itemIm);
	}

	if (n != 0) {
		ui.tableWidget->selectRow(0);
		ui.buttonCopyImage->setEnabled(true);
		ui.buttonCopyXml->setEnabled(true);
	}
	else {
		clearPreview();
		ui.buttonCopyImage->setEnabled(false);
		ui.buttonCopyXml->setEnabled(false);
	}

	QApplication::restoreOverrideCursor();
}

static QString toParamString(const QString& name, const complex& z) {
	return "Re_" + name + "=\"" + QString::number(z.real(), 'g', 16) + "\" "
		"Im_" + name + "=\"" + QString::number(z.imag(), 'g', 16) + "\" ";
}

static QString toXformString(const QString& name, const Moebius<complex>& M) {
	return "<xform name=\"" + name + "\" "
		"weight=\"1\" color=\"0\" mobius=\"1\" coefs=\"1 0 0 1 0 0\" " +
		toParamString("A", M.a) +
		toParamString("B", M.b) +
		toParamString("C", M.c) +
		toParamString("D", M.d) +
		"/>\n";
}

void MainForm::copyImageButtonClicked() {
	QApplication::clipboard()->setImage(previewImage);
}

void MainForm::copyXmlButtonClicked() {
	QString viewScale, viewXform;
	int calcMode = ui.comboBoxCalcMode->currentIndex();
	int viewMode = ui.comboBoxViewMode->currentIndex();
	viewScale = (viewMode == 0) ? "512" : "256";
	if (viewMode != calcMode) {
		static const QString views[] = {
			"Re_A=\"1\" Im_A=\"0\" Re_B=\"-1\" Im_B=\"0\" Re_C=\"1\" Im_C=\"0\" Re_D=\"1\" Im_D=\"0\"",
			"Re_A=\"1\" Im_A=\"0\" Re_B=\"1\" Im_B=\"0\" Re_C=\"-1\" Im_C=\"0\" Re_D=\"1\" Im_D=\"0\""
		};
		viewXform = "<finalxform color=\"0\" symmetry=\"1\" mobius=\"1\" coefs=\"1 0 0 1 0 0\" " +
					views[viewMode] + "/>\n";
	}
	auto nStr = [](int n) { return isInf(n) ? "inf" : QString::number(n); };
	QString kleinianParams =
		QString::number(K.a) + " " + QString::number(K.b) + " " +
		nStr(K.n) + " " + nStr(K.h) + " " + nStr(K.v1) + " " + nStr(K.v2);

	QString xml =
		"<flame name=\"Kleinian " + kleinianParams + "\""
		" size=\"1024 1024\" center=\"0 0\" scale=\"" + viewScale + "\""
		" oversample=\"1\" filter=\"0.5\" quality=\"50\""
		" background=\"0 0 0\" brightness=\"4\" gamma=\"4\" gamma_threshold=\"0.04\">\n";

	if (ui.groupBox_a->isChecked())
		xml += toXformString("a", Ma);
	if (ui.groupBox_A->isChecked())
		xml += toXformString("A", Ma.inverse());
	if (ui.groupBox_b->isChecked())
		xml += toXformString("b", Mb);
	if (ui.groupBox_B->isChecked())
		xml += toXformString("B", Mb.inverse());
	if (ui.groupBox_c->isChecked())
		xml += toXformString("c", Mc);
	if (ui.groupBox_C->isChecked())
		xml += toXformString("C", Mc.inverse());
	if (ui.groupBox_x->isChecked())
		xml += toXformString(ui.lineEditSeq1->text(), Mseq1);
	if (ui.groupBox_X->isChecked())
		xml += toXformString(ui.lineEditSeq2->text(), Mseq2);

	xml += viewXform +
		"<color index=\"0\" rgb=\"255 255 255\"/>\n"
		"</flame>";

	QApplication::clipboard()->setText(xml);
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

void MainForm::valueHbuttonClicked() {
	QApplication::clipboard()->setText(QString::fromStdString(boost::lexical_cast<std::string>(offsetH)));
}

void MainForm::valueV1buttonClicked() {
	QApplication::clipboard()->setText(QString::fromStdString(boost::lexical_cast<std::string>(offsetV1)));
}

void MainForm::valueV2buttonClicked() {
	QApplication::clipboard()->setText(QString::fromStdString(boost::lexical_cast<std::string>(offsetV2)));
}

static void displayMatrix(const Moebius<complex>& M, QLabel* la, QLabel* lb, QLabel* lc, QLabel* ld) {
	auto toString = [](complex z){
		return (z.imag() != 0)
			? QString::number(z.real(), 'g', 8) + "\n" + QString::number(z.imag(), 'g', 8)
			: QString::number(z.real(), 'g', 8);
	};
	la->setText(toString(M.a));
	lb->setText(toString(M.b));
	lc->setText(toString(M.c));
	ld->setText(toString(M.d));
}

void MainForm::tableItemSelectionChanged() {
	auto selection = ui.tableWidget->selectedItems();
	if (selection.isEmpty()) {
		return;
	}

	int row = selection.front()->row();
	if (row >= 0) {
		auto itemRe = ui.tableWidget->item(row, 0);
		auto itemIm = ui.tableWidget->item(row, 1);
		if (itemRe && itemIm) {
			double re = itemRe->data(Qt::UserRole).value<double>();
			double im = itemIm->data(Qt::UserRole).value<double>();

			auto calcMode = qBound(0, ui.comboBoxCalcMode->currentIndex(), 1);
			auto transforms = K.createGenerators(complex(re, im), calcMode == 0);
			Ma = transforms[0];
			Mb = transforms[1];
			Mc = transforms[2];

			displayMatrix(Ma, ui.label_a11, ui.label_a12, ui.label_a21, ui.label_a22);
			displayMatrix(Mb, ui.label_b11, ui.label_b12, ui.label_b21, ui.label_b22);
			displayMatrix(Mc, ui.label_c11, ui.label_c12, ui.label_c21, ui.label_c22);
			displayMatrix(Ma.inverse(), ui.label_A11, ui.label_A12, ui.label_A21, ui.label_A22);
			displayMatrix(Mb.inverse(), ui.label_B11, ui.label_B12, ui.label_B21, ui.label_B22);
			displayMatrix(Mc.inverse(), ui.label_C11, ui.label_C12, ui.label_C21, ui.label_C22);

			updateMseq1(ui.lineEditSeq1->text());
			updateMseq2(ui.lineEditSeq2->text());

			renderPreview();
			return;
		}
	}

	clearPreview();
}

void MainForm::updateMseq1(const QString& str) {
	Mseq1 = sequenceMatrix(str.toStdString());
	displayMatrix(Mseq1, ui.label_x11, ui.label_x12, ui.label_x21, ui.label_x22);
}

void MainForm::updateMseq2(const QString& str) {
	Mseq2 = sequenceMatrix(str.toStdString());
	displayMatrix(Mseq2, ui.label_X11, ui.label_X12, ui.label_X21, ui.label_X22);
}

void MainForm::lineEditSeq1Edited(const QString& str) {
	updateMseq1(str);
	renderPreview();
}

void MainForm::lineEditSeq2Edited(const QString& str) {
	updateMseq2(str);
	renderPreview();
}

void MainForm::clearPreview() {
	QPixmap pixmap(renderWidth, renderHeight);
	pixmap.fill(Qt::black);
	ui.labelPreview->setPixmap(pixmap);
}

void MainForm::renderPreview() {
	std::vector<Moebius<Renderer::complex>> transforms;

	if (ui.groupBox_a->isChecked())
		transforms.push_back(Ma);
	if (ui.groupBox_A->isChecked())
		transforms.push_back(Ma.inverse());
	if (ui.groupBox_b->isChecked())
		transforms.push_back(Mb);
	if (ui.groupBox_B->isChecked())
		transforms.push_back(Mb.inverse());
	if (ui.groupBox_c->isChecked())
		transforms.push_back(Mc);
	if (ui.groupBox_C->isChecked())
		transforms.push_back(Mc.inverse());
	if (ui.groupBox_x->isChecked())
		transforms.push_back(Mseq1);
	if (ui.groupBox_X->isChecked())
		transforms.push_back(Mseq2);

	auto calcMode = (Renderer::ViewMode)qBound(0, ui.comboBoxCalcMode->currentIndex(), 1);
	auto viewMode = (Renderer::ViewMode)qBound(0, ui.comboBoxViewMode->currentIndex(), 1);

	renderer.clear();
	if (!transforms.empty())
		renderer.iterate(renderIterations, transforms, calcMode, viewMode);
	renderer.tonemap(previewImage);
	ui.labelPreview->setPixmap(QPixmap::fromImage(previewImage));
}
