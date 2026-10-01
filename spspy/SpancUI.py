from .Spanc import Spanc, PeakType, DEG2RAD
from .SpancExport import export_energy_calibration_csv, export_spanc_workbook
from .Fitter import convert_fit_points_to_arrays, convert_resid_points_to_arrays
from .ui.MPLCanvas import MPLCanvas
from .ui.ReactionDialog import ReactionDialog
from .ui.TargetDialog import TargetDialog
from .ui.PeakDialog import PeakDialog

from PySide6.QtWidgets import QApplication, QWidget, QMainWindow
from PySide6.QtWidgets import QLabel, QTabWidget, QTableWidget, QTableWidgetItem
from PySide6.QtWidgets import QHBoxLayout, QVBoxLayout, QGroupBox
from PySide6.QtWidgets import QPushButton, QTextEdit, QSpinBox, QComboBox
from PySide6.QtWidgets import QDoubleSpinBox, QFileDialog, QMessageBox
from PySide6.QtGui import QAction

# from qdarktheme import load_stylesheet, load_palette
import matplotlib as mpl
import numpy as np
from numpy.typing import NDArray
import sys
import pickle

#Get y-value for baseline at 0
def baseline(x: float) -> float:
    return 0.0

class SpancGUI(QMainWindow):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.setWindowTitle("SPANC")
        self.spanc = Spanc()
        self.tablelayout = QVBoxLayout()
        self.plotlayout = QHBoxLayout()
        self.layout = QVBoxLayout()
        self.centralWidget = QTabWidget(self)
        self.setCentralWidget(self.centralWidget)
        self.centralWidget.setLayout(self.layout)
        self.tableTab = QWidget(self.centralWidget)
        self.tableTab.setLayout(self.tablelayout)
        self.plotTab = QWidget(self.centralWidget)
        self.plotTab.setLayout(self.plotlayout)
        self.energyTab = QWidget(self.centralWidget)
        self.energyLayout = QHBoxLayout()
        self.energyTab.setLayout(self.energyLayout)
        self.centralWidget.addTab(self.tableTab, "Data Tables")
        self.centralWidget.addTab(self.plotTab, "Plots and Fits")
        self.centralWidget.addTab(self.energyTab, "Energy Calibration")
        self.fitFlag = False
        self.create_menus()
        self.create_fit_canvas()
        self.create_target_table() 
        self.create_reaction_table()
        self.create_calibration_table()
        self.create_output_table()
        self.create_fit_result_text()
        self.create_energy_calibration_tab()
        self.refresh_energy_reaction_box()
        self.show()

    def create_menus(self) -> None:
        self.fileMenu = self.menuBar().addMenu("&File")
        saveAction = QAction("&Save...",self)
        openAction = QAction("&Open...",self)
        saveFitAction = QAction("Save Calibration Plot...", self)
        saveResidualAction = QAction("Save Residual Plot...", self)
        exportWorkbookAction = QAction("Export SPANC Workbook...", self)
        self.fileMenu.addAction(saveAction)
        self.fileMenu.addAction(openAction)
        self.fileMenu.addAction(saveFitAction)
        self.fileMenu.addAction(saveResidualAction)
        self.fileMenu.addAction(exportWorkbookAction)
        self.fileMenu.addAction("&Exit", self.close)
        saveAction.triggered.connect(self.handle_save)
        openAction.triggered.connect(self.handle_open)
        saveFitAction.triggered.connect(self.handle_save_fit)
        saveResidualAction.triggered.connect(self.handle_save_residual)
        exportWorkbookAction.triggered.connect(self.handle_export_workbook)
        
        self.addMenu = self.menuBar().addMenu("&New")
        newTargetAction = QAction("New target...", self)
        newReactionAction = QAction("New reaction...", self)
        newCalibrationAction = QAction("New calibration...", self)
        newOutputAction = QAction("New output...", self)
        self.addMenu.addAction(newTargetAction)
        self.addMenu.addAction(newReactionAction)
        self.addMenu.addAction(newCalibrationAction)
        self.addMenu.addAction(newOutputAction)
        newTargetAction.triggered.connect(self.handle_new_target)
        newReactionAction.triggered.connect(self.handle_new_reaction)
        newCalibrationAction.triggered.connect(self.handle_new_calibration)
        newOutputAction.triggered.connect(self.handle_new_output)

    def create_fit_canvas(self) -> None:
        self.fitGroup = QGroupBox("Calibration Fit", self.plotTab)
        fitLayout = QVBoxLayout()
        self.fitCanvas = MPLCanvas(self.fitGroup, width=6, height=6, dpi=100)
        self.residCanvas = MPLCanvas(self.fitGroup, width=6, height=6, dpi=100)

        self.fitOptionGroup = QGroupBox("Fit options", self.fitGroup)
        fitOptionLayout = QHBoxLayout()
        self.fitButton = QPushButton("Run Fit", self.fitOptionGroup)
        self.fitButton.clicked.connect(self.handle_run_fit)
        self.fitOrderBox = QSpinBox(self.fitOptionGroup)
        self.fitOrderBox.valueChanged.connect(self.handle_change_fit_order)
        self.fitOrderBox.setValue(1)
        self.fitOrderBox.setMaximum(10)
        self.fitOrderBox.setMinimum(0)
        fitOptionLayout.addWidget(QLabel("Polynomial Order", self.fitOptionGroup))
        fitOptionLayout.addWidget(self.fitOrderBox)
        fitOptionLayout.addWidget(self.fitButton)
        self.fitOptionGroup.setLayout(fitOptionLayout)

        fitLayout.addWidget(QLabel("Fit", self.fitCanvas))
        fitLayout.addWidget(self.fitCanvas)
        fitLayout.addWidget(QLabel("Residuals", self.fitCanvas))
        fitLayout.addWidget(self.residCanvas)
        fitLayout.addWidget(self.fitOptionGroup)
        self.fitGroup.setLayout(fitLayout)
        self.plotlayout.addWidget(self.fitGroup)

    def create_target_table(self) -> None:
        self.targetGroup = QGroupBox("Targets", self.tableTab)
        targetLayout = QVBoxLayout()
        self.targetTable = QTableWidget(self.targetGroup)
        self.targetTable.setColumnCount(6)
        self.targetTable.setHorizontalHeaderLabels(["L1 Thickness(ug/cm^2)", "L1 Compound","L2 Thickness(ug/cm^2)", "L2 Compound","L3 Thickness(ug/cm^2)", "L3 Compound"])
        targetLayout.addWidget(self.targetTable)
        self.targetGroup.setLayout(targetLayout)
        self.tablelayout.addWidget(self.targetGroup)
        self.targetTable.resizeColumnsToContents()
        self.targetTable.cellDoubleClicked.connect(self.handle_update_target)
        self.targetTable.setEditTriggers(QTableWidget.EditTrigger.NoEditTriggers)

    def create_reaction_table(self) -> None:
        self.rxnGroup = QGroupBox("Reactions", self.tableTab)
        rxnLayout = QVBoxLayout()
        self.reactionTable = QTableWidget(self.rxnGroup)
        self.reactionTable.setColumnCount(5)
        self.reactionTable.setHorizontalHeaderLabels(["Target Material","Reaction Equation","Beam KE(MeV)","BField(kG)","Angle(deg)"])
        rxnLayout.addWidget(self.reactionTable)
        self.rxnGroup.setLayout(rxnLayout)
        self.tablelayout.addWidget(self.rxnGroup)
        self.reactionTable.resizeColumnsToContents()
        self.reactionTable.cellDoubleClicked.connect(self.handle_update_reaction)
        self.reactionTable.setEditTriggers(QTableWidget.EditTrigger.NoEditTriggers)

    def create_calibration_table(self) -> None:
        self.calGroup = QGroupBox("Calibration Peaks", self.tableTab)
        calLayout = QVBoxLayout()
        self.calibrationTable = QTableWidget(self.calGroup)
        self.calibrationTable.setColumnCount(9)
        self.calibrationTable.setHorizontalHeaderLabels(["Peak ID","Reaction","x(mm)","ux stat.(mm)","ux sys.(mm)","rho(cm)","urho(cm)","Ex(MeV)","uEx(MeV)"])
        calLayout.addWidget(self.calibrationTable)
        self.calGroup.setLayout(calLayout)
        self.tablelayout.addWidget(self.calGroup)
        self.calibrationTable.resizeColumnsToContents()
        self.calibrationTable.cellDoubleClicked.connect(self.handle_update_calibration)
        self.calibrationTable.setEditTriggers(QTableWidget.EditTrigger.NoEditTriggers)

    def create_output_table(self) -> None:
        self.outGroup = QGroupBox("Output Peaks", self.tableTab)
        outLayout = QVBoxLayout()
        self.outputTable = QTableWidget(self.outGroup)
        self.outputTable.setColumnCount(13)
        self.outputTable.setHorizontalHeaderLabels(["Peak ID", "Reaction","x(mm)","ux stat.(mm)","ux sys.(mm)","rho(cm)","urho(cm)","Ex(MeV)","uEx(MeV)","FWHM(mm)","uFWHM(mm)","FWHM(MeV)","uFWHM(MeV)"])
        outLayout.addWidget(self.outputTable)
        self.outGroup.setLayout(outLayout)
        self.tablelayout.addWidget(self.outGroup)
        self.outputTable.resizeColumnsToContents()
        self.outputTable.cellDoubleClicked.connect(self.handle_update_output)
        self.outputTable.setEditTriggers(QTableWidget.EditTrigger.NoEditTriggers)

    def create_fit_result_text(self) -> None:
        self.fitTextGroup = QGroupBox("Fit Results", self.plotTab)
        fitTextLayout = QVBoxLayout()
        self.fitResultText = QTextEdit(self.fitTextGroup)
        self.fitResultText.setReadOnly(True)
        fitTextLayout.addWidget(self.fitResultText)
        self.fitTextGroup.setLayout(fitTextLayout)
        self.plotlayout.addWidget(self.fitTextGroup)

    def create_energy_calibration_tab(self) -> None:
        self.energyPlotGroup = QGroupBox("Direct Position-to-Energy Fit", self.energyTab)
        energyPlotLayout = QVBoxLayout()

        self.energyCanvas = MPLCanvas(
            self.energyPlotGroup, width=7, height=6, dpi=100
        )
        energyPlotLayout.addWidget(self.energyCanvas)

        self.energyOptionGroup = QGroupBox("Energy fit options", self.energyPlotGroup)
        energyOptionLayout = QHBoxLayout()
        self.energyReactionBox = QComboBox(self.energyOptionGroup)
        self.energyReactionBox.currentIndexChanged.connect(
            self.handle_energy_reaction_changed
        )
        self.energyFitOrderBox = QSpinBox(self.energyOptionGroup)
        self.energyFitOrderBox.setRange(0, 10)
        self.energyFitOrderBox.setValue(1)
        self.energyFitOrderBox.valueChanged.connect(
            self.handle_change_energy_fit_order
        )
        self.energyFitButton = QPushButton("Run Energy Fit", self.energyOptionGroup)
        self.energyFitButton.clicked.connect(self.handle_run_energy_fit)
        self.saveEnergyPlotButton = QPushButton(
            "Save Plot...", self.energyOptionGroup
        )
        self.saveEnergyPlotButton.clicked.connect(self.handle_save_energy_plot)
        self.exportEnergyButton = QPushButton(
            "Export CSV...", self.energyOptionGroup
        )
        self.exportEnergyButton.clicked.connect(
            self.handle_export_energy_calibration
        )

        energyOptionLayout.addWidget(QLabel("Reaction", self.energyOptionGroup))
        energyOptionLayout.addWidget(self.energyReactionBox)
        energyOptionLayout.addWidget(
            QLabel("Polynomial Order", self.energyOptionGroup)
        )
        energyOptionLayout.addWidget(self.energyFitOrderBox)
        energyOptionLayout.addWidget(self.energyFitButton)
        energyOptionLayout.addWidget(self.saveEnergyPlotButton)
        energyOptionLayout.addWidget(self.exportEnergyButton)
        self.energyOptionGroup.setLayout(energyOptionLayout)
        energyPlotLayout.addWidget(self.energyOptionGroup)
        self.energyPlotGroup.setLayout(energyPlotLayout)
        self.energyLayout.addWidget(self.energyPlotGroup)

        self.energyResultGroup = QGroupBox("Energy Fit Results", self.energyTab)
        energyResultLayout = QVBoxLayout()
        self.energyResultText = QTextEdit(self.energyResultGroup)
        self.energyResultText.setReadOnly(True)
        energyResultLayout.addWidget(self.energyResultText)

        calculatorGroup = QGroupBox("Position-to-Energy Calculator", self.energyResultGroup)
        calculatorLayout = QVBoxLayout()
        calculatorInputLayout = QHBoxLayout()
        self.energyPositionBox = QDoubleSpinBox(calculatorGroup)
        self.energyPositionBox.setRange(-1000000.0, 1000000.0)
        self.energyPositionBox.setDecimals(6)
        self.energyPositionBox.setSingleStep(0.1)
        self.energyPositionBox.valueChanged.connect(self.update_energy_calculator)
        calculatorInputLayout.addWidget(QLabel("Position (mm)", calculatorGroup))
        calculatorInputLayout.addWidget(self.energyPositionBox)
        calculatorLayout.addLayout(calculatorInputLayout)
        self.energyCalculatorResult = QLabel(
            "Run an energy fit to enable the calculator.", calculatorGroup
        )
        self.energyCalculatorResult.setWordWrap(True)
        calculatorLayout.addWidget(self.energyCalculatorResult)
        calculatorGroup.setLayout(calculatorLayout)
        energyResultLayout.addWidget(calculatorGroup)

        self.energyResultGroup.setLayout(energyResultLayout)
        self.energyLayout.addWidget(self.energyResultGroup)

    def handle_save(self) -> None:
        fileName = QFileDialog.getSaveFileName(self, "Save Input","./","SPANC Files (*.spanc)")
        if fileName[0]:
            #self.spanc.WriteConfig(fileName[0])
            with open(fileName[0], "wb") as savefile:
                pickle.dump(self.spanc, savefile, pickle.HIGHEST_PROTOCOL)
                savefile.close()

    def handle_save_fit(self) -> None:
        fileName = QFileDialog.getSaveFileName(
            self, "Save Calibration Plot", "./", "Image Files (*.png *.eps *.pdf)"
        )
        if fileName[0]:
            self.fitCanvas.fig.savefig(fileName[0])

    def handle_save_residual(self) -> None:
        fileName = QFileDialog.getSaveFileName(
            self, "Save Residual Plot", "./", "Image Files (*.png *.eps *.pdf)"
        )
        if fileName[0]:
            self.residCanvas.fig.savefig(fileName[0])

    def handle_export_workbook(self) -> None:
        fileName = QFileDialog.getSaveFileName(
            self,
            "Export SPANC Workbook",
            "./SPANC_analysis.xlsx",
            "Excel Workbooks (*.xlsx)",
        )
        if not fileName[0]:
            return
        try:
            savedPath = export_spanc_workbook(self.spanc, fileName[0])
            self.update_output_table()
            QMessageBox.information(
                self, "SPANC Export", f"Workbook saved to:\n{savedPath}"
            )
        except Exception as error:
            QMessageBox.warning(
                self, "SPANC Export Failed", f"Could not export workbook:\n{error}"
            )

    def handle_save_energy_plot(self) -> None:
        rxnName = self.get_selected_energy_reaction()
        if rxnName is None or not self.spanc.has_energy_fit(rxnName):
            QMessageBox.warning(
                self,
                "No Energy Fit",
                "Run an energy calibration fit before saving its plot.",
            )
            return
        defaultName = f"{self.energy_export_stem(rxnName)}.png"
        fileName = QFileDialog.getSaveFileName(
            self,
            "Save Energy Calibration Plot",
            f"./{defaultName}",
            "Image Files (*.png *.eps *.pdf)",
        )
        if fileName[0]:
            self.energyCanvas.fig.savefig(fileName[0])

    def handle_export_energy_calibration(self) -> None:
        rxnName = self.get_selected_energy_reaction()
        if rxnName is None or not self.spanc.has_energy_fit(rxnName):
            QMessageBox.warning(
                self,
                "No Energy Fit",
                "Run an energy calibration fit for the selected reaction "
                "before exporting it.",
            )
            return

        defaultName = f"{self.energy_export_stem(rxnName)}.csv"
        fileName = QFileDialog.getSaveFileName(
            self,
            "Export Selected Energy Calibration",
            f"./{defaultName}",
            "CSV Files (*.csv)",
        )
        if not fileName[0]:
            return
        try:
            savedPath = export_energy_calibration_csv(
                self.spanc, rxnName, fileName[0]
            )
            QMessageBox.information(
                self, "Energy Calibration Export", f"CSV saved to:\n{savedPath}"
            )
        except Exception as error:
            QMessageBox.warning(
                self,
                "Energy Calibration Export Failed",
                f"Could not export energy calibration:\n{error}",
            )

    def handle_open(self) -> None:
        fileName = QFileDialog.getOpenFileName(self, "Open Input","./","SPANC Files (*.spanc)")
        if fileName[0]:
            with open(fileName[0], "rb") as openfile:
                self.spanc = pickle.load(openfile)
                self.spanc.ensure_energy_fitters()
                self.update_target_table()
                self.update_reaction_table()
                self.update_calibration_table()
                self.update_fit_order()
                self.update_output_table()
                self.refresh_energy_reaction_box()
                openfile.close()

    def handle_new_target(self) -> None:
        targDia = TargetDialog(self)
        targDia.new_target.connect(self.spanc.add_target)
        if targDia.exec() :
            self.update_target_table()
        return

    def handle_update_target(self, row: int, col: int) -> None:
        targName = self.targetTable.verticalHeaderItem(row).text()
        targDia = TargetDialog(self, target=self.spanc.targets[targName])
        targDia.new_target.connect(self.spanc.add_target)
        if targDia.exec():
            self.update_target_table()
            self.spanc.calculate_calibrations()
            self.update_reaction_table()
            self.update_calibration_table()
            self.spanc.calculate_outputs()
            self.update_output_table()
        return

    def handle_new_reaction(self) -> None:
        rxnDia = ReactionDialog(self, targets=self.spanc.targets.keys(), extraParams=True)
        rxnDia.new_reaction.connect(self.spanc.add_reaction)
        if rxnDia.exec():
            self.update_reaction_table()
            self.refresh_energy_reaction_box()
        return

    def handle_update_reaction(self, row: int, col: int) -> None:
        rxnName = self.reactionTable.verticalHeaderItem(row).text()
        rxnDia = ReactionDialog(self, targets=self.spanc.targets.keys(), rxn=self.spanc.reactions[rxnName], rxnKey=rxnName, extraParams=True)
        rxnDia.update_reaction.connect(self.spanc.update_reaction_parameters)
        if rxnDia.exec():
            self.update_reaction_table()
            self.refresh_energy_reaction_box()
            self.spanc.calculate_calibrations()
            self.update_calibration_table()
            self.spanc.calculate_outputs()
            self.update_output_table()
        return

    def handle_new_calibration(self) -> None:
        calDia = PeakDialog(PeakType.CALIBRATION, self.spanc.reactions.keys(), self)
        calDia.new_peak.connect(self.spanc.add_calibration)
        if calDia.exec():
            self.update_calibration_table()
            self.refresh_energy_display()
        return
        
    def handle_update_calibration(self, row: int, col: int) -> None:
        peakID = int(self.calibrationTable.item(row, 0).text())
        peakData = self.spanc.calibrations[peakID]
        calDia = PeakDialog(PeakType.CALIBRATION, self.spanc.reactions.keys(), self, peak=peakData)
        calDia.new_peak.connect(self.spanc.add_calibration)
        calDia.delete_peak.connect(self.spanc.remove_calibration)
        if calDia.exec():
            self.update_calibration_table()
            self.refresh_energy_display()
            if self.spanc.isFit == True:
                self.handle_run_fit()
        return

    def handle_new_output(self) -> None:
        outDia = PeakDialog(PeakType.OUTPUT, self.spanc.reactions.keys(), self)
        outDia.new_peak.connect(self.spanc.add_output)
        if outDia.exec():
            self.update_output_table()
        return

    def handle_update_output(self, row: int, col: int) -> None:
        peakID = int(self.calibrationTable.item(row, 0).text())
        peakData = self.spanc.outputs[peakID]
        outDia = PeakDialog(PeakType.OUTPUT, self.spanc.reactions.keys(), self, peak=peakData)
        outDia.new_peak.connect(self.spanc.add_output)
        if outDia.exec():
            self.update_output_table()
        return

    def handle_change_fit_order(self, order: int) -> None:
        self.spanc.set_fit_order(order)        

    def get_selected_energy_reaction(self) -> str | None:
        return self.energyReactionBox.currentData()

    def energy_export_stem(self, rxnName: str) -> str:
        angle = self.spanc.reactions[rxnName].params.spsAngle / DEG2RAD
        angleText = f"{angle:.6g}".replace("-", "m").replace(".", "p")
        return f"{rxnName}_{angleText}deg_energy_calibration"

    def refresh_energy_reaction_box(self) -> None:
        previousReaction = self.get_selected_energy_reaction()
        self.energyReactionBox.blockSignals(True)
        self.energyReactionBox.clear()
        selectedIndex = -1
        for index, (rxnName, rxn) in enumerate(self.spanc.reactions.items()):
            angle = rxn.params.spsAngle / DEG2RAD
            label = (
                f"{rxnName} - {rxn}, angle={angle:.6g} deg, "
                f"B={rxn.params.magneticField:.6g} kG"
            )
            self.energyReactionBox.addItem(label, rxnName)
            if rxnName == previousReaction:
                selectedIndex = index
        if self.energyReactionBox.count() > 0:
            self.energyReactionBox.setCurrentIndex(
                selectedIndex if selectedIndex >= 0 else 0
            )
        self.energyReactionBox.blockSignals(False)
        self.handle_energy_reaction_changed()

    def refresh_energy_display(self) -> None:
        self.handle_energy_reaction_changed()

    def clear_energy_display(self, message: str) -> None:
        self.energyCanvas.axes.cla()
        self.energyCanvas.axes.set_xlabel("Focal-Plane Position (mm)")
        self.energyCanvas.axes.set_ylabel("Excitation Energy (MeV)")
        self.energyCanvas.fig.tight_layout()
        self.energyCanvas.draw()
        self.energyResultText.setMarkdown(f"# Energy Calibration\n\n{message}")
        self.energyCalculatorResult.setText(
            "Run an energy fit to enable the calculator."
        )

    def handle_energy_reaction_changed(self, index: int = -1) -> None:
        rxnName = self.get_selected_energy_reaction()
        if rxnName is None:
            self.clear_energy_display(
                "Add a reaction before running a direct energy calibration."
            )
            return

        fitter = self.spanc.get_energy_fitter(rxnName)
        self.energyFitOrderBox.blockSignals(True)
        self.energyFitOrderBox.setValue(fitter.polynomialOrder)
        self.energyFitOrderBox.blockSignals(False)

        if self.spanc.has_energy_fit(rxnName):
            self.draw_energy_fit(rxnName)
            self.update_energy_fit_text(rxnName)
            self.update_energy_calculator()
        else:
            pointCount = len(self.spanc.get_calibrations_for_reaction(rxnName))
            self.clear_energy_display(
                f"Selected **{rxnName}** with {pointCount} calibration point(s). "
                "Choose the polynomial order and run the energy fit."
            )

    def handle_change_energy_fit_order(self, order: int) -> None:
        rxnName = self.get_selected_energy_reaction()
        if rxnName is None:
            return
        self.spanc.set_energy_fit_order(rxnName, order)
        self.clear_energy_display(
            f"The polynomial order for **{rxnName}** is now {order}. "
            "Run the energy fit to update the calibration."
        )

    def handle_run_energy_fit(self) -> None:
        rxnName = self.get_selected_energy_reaction()
        if rxnName is None:
            QMessageBox.warning(
                self, "No Reaction", "Add and select a reaction before fitting."
            )
            return
        try:
            self.spanc.fit_energy(rxnName)
        except ValueError as error:
            QMessageBox.warning(self, "Energy Fit Cannot Run", str(error))
            return

        self.draw_energy_fit(rxnName)
        self.update_energy_fit_text(rxnName)
        self.update_energy_calculator()

        warnings = self.spanc.get_energy_fitter(rxnName).lastWarnings
        if warnings:
            QMessageBox.warning(self, "Energy Fit Uncertainty Warning", "\n\n".join(warnings))

    def draw_energy_fit(self, rxnName: str) -> None:
        fitter = self.spanc.get_energy_fitter(rxnName)
        fitData = fitter.fitData or []
        if not fitData:
            self.clear_energy_display("No energy-fit data are available.")
            return

        xArray, yArray, xErrArray, yErrArray = convert_fit_points_to_arrays(
            fitData
        )
        xMin = float(np.amin(xArray))
        xMax = float(np.amax(xArray))
        padding = max((xMax - xMin) * 0.05, 1.0)
        fitArray = np.linspace(xMin - padding, xMax + padding, 1000)

        self.energyCanvas.axes.cla()
        self.energyCanvas.axes.errorbar(
            xArray,
            yArray,
            xerr=xErrArray,
            yerr=yErrArray,
            marker="o",
            linestyle="None",
            elinewidth=2.0,
            capsize=3.0,
            label="Calibration peaks",
        )
        self.energyCanvas.axes.plot(
            fitArray,
            fitter.evaluate(fitArray),
            label=f"Order {fitter.polynomialOrder} fit",
        )
        self.energyCanvas.axes.set_xlabel("Focal-Plane Position (mm)")
        self.energyCanvas.axes.set_ylabel("Excitation Energy (MeV)")
        self.energyCanvas.axes.legend()
        self.energyCanvas.fig.tight_layout()
        self.energyCanvas.draw()

    def update_energy_fit_text(self, rxnName: str) -> None:
        fitter = self.spanc.get_energy_fitter(rxnName)
        rxn = self.spanc.reactions[rxnName]
        angle = rxn.params.spsAngle / DEG2RAD
        equationMeV, equationKeV = self.spanc.get_energy_equations(rxnName)
        residualsKeV = self.spanc.get_direct_energy_residuals(rxnName)

        parameterLines = []
        for order, (value, uncertainty) in enumerate(
            zip(fitter.get_parameters(), fitter.get_parameter_errors())
        ):
            if order == 0:
                units = "MeV"
            elif order == 1:
                units = "MeV/mm"
            else:
                units = f"MeV/mm^{order}"
            parameterLines.append(
                f"- b{order} = {value:.10g} +/- {uncertainty:.6g} {units}"
            )

        warningText = ""
        if fitter.lastWarnings:
            warningText = "## Fit Warnings\n\n" + "\n\n".join(
                f"- {warning}" for warning in fitter.lastWarnings
            ) + "\n\n"

        markdownString = (
            f"# Energy Calibration: {rxnName}\n\n"
            f"**Reaction:** {rxn}\n\n"
            f"**Angle:** {angle:.8g} deg  \n"
            f"**Magnetic Field:** {rxn.params.magneticField:.8g} kG  \n"
            f"**Beam Energy:** {rxn.params.beamEnergy:.8g} MeV\n\n"
            f"## Copyable Equations\n\n"
            f"```text\n{equationMeV}\n{equationKeV}\n```\n\n"
            f"## Fit Statistics\n\n"
            f"- Polynomial order: {fitter.polynomialOrder}\n"
            f"- Calibration points: {len(fitter.fitData or [])}\n"
            f"- Chi-square: {fitter.get_chisquare():.8g}\n"
            f"- NDF: {fitter.get_ndf()}\n"
            f"- Reduced chi-square: {fitter.get_reduced_chisquare():.8g}\n"
            f"- R-squared (vertical-residual diagnostic): "
            f"{fitter.get_r_squared():.10g}\n\n"
            f"## Parameters\n\n"
            + "\n".join(parameterLines)
            + "\n\n"
            f"## Direct Energy Residuals (keV)\n\n"
            f"{np.array_str(residualsKeV, precision=4)}\n\n"
            + warningText
        )
        self.energyResultText.setMarkdown(markdownString)

    def update_energy_calculator(self, position: float | None = None) -> None:
        rxnName = self.get_selected_energy_reaction()
        if rxnName is None or not self.spanc.has_energy_fit(rxnName):
            self.energyCalculatorResult.setText(
                "Run an energy fit to enable the calculator."
            )
            return
        if position is None:
            position = self.energyPositionBox.value()

        excitationMeV = self.spanc.evaluate_energy(rxnName, position)
        sensitivity = self.spanc.get_energy_sensitivity(rxnName, position)
        self.energyCalculatorResult.setText(
            f"Ex = {excitationMeV:.10g} MeV = "
            f"{excitationMeV * 1000.0:.10g} keV\n"
            f"dEx/dx = {sensitivity:.10g} keV/mm"
        )

    def handle_run_fit(self) -> None:
        order = self.spanc.fitter.polynomialOrder
        npoints = len(self.spanc.calibrations)
        if npoints < (order + 2):
            print(f"Warning! Attempting to fit {npoints} data points with order {order} polyomial, too few degrees of freedom!")
            print(f"Increase number of data points to at minimum {order+2} to use a polynomial of this order.")
            return
        fitData = self.spanc.fit()
        xArray, yArray, xErrArray, yErrArray = convert_fit_points_to_arrays(fitData)
        xMin = np.amin(xArray)
        xMax = np.amax(xArray)
        fitArray = np.linspace(xMin, xMax, 1000)
        self.fitCanvas.axes.cla()
        self.fitCanvas.axes.errorbar(xArray, yArray, yerr=yErrArray, xerr=xErrArray, marker="o", linestyle="None", elinewidth=2.0)
        self.fitCanvas.axes.plot(fitArray, self.spanc.fitter.evaluate(fitArray))
        self.fitCanvas.axes.set_xlabel(r"$x$ (mm)")
        self.fitCanvas.axes.set_ylabel(r"$\rho$ (cm)")
        self.fitCanvas.fig.tight_layout()
        self.fitCanvas.draw()
        self.spanc.calculate_outputs()
        self.update_output_table()
        self.fitFlag = True

        residData = self.spanc.get_residuals()
        xArray, residArray, residErrArray, studentResidArray = convert_resid_points_to_arrays(residData)
        energyResidArray, energyResidErrArray = self.spanc.get_energy_residuals()
        self.residCanvas.axes.cla()
        self.residCanvas.axes.errorbar(
            xArray,
            residArray,
            yerr=residErrArray,
            marker="o",
            linestyle="None",
            elinewidth=2.0,
            capsize=3.0,
        )
        self.residCanvas.axes.hlines(0.0, xMin, xMax, colors="r", linestyles="dashed")
        self.residCanvas.axes.set_xlabel(r"$x$ (mm)")
        self.residCanvas.axes.set_ylabel(r"Residual (cm)")
        self.residCanvas.fig.tight_layout()
        self.residCanvas.draw()
        self.update_fit_text(
            residArray,
            studentResidArray,
            energyResidArray,
            energyResidErrArray,
        )

    def update_target_table(self) -> None:
        self.targetTable.setRowCount(len(self.spanc.targets))
        self.targetTable.setVerticalHeaderLabels(self.spanc.targets.keys())
        for row, key in enumerate(self.spanc.targets):
            for i, layer in enumerate(self.spanc.targets[key].layer_details) :
                self.targetTable.setItem(row, 0+i*2, QTableWidgetItem(str(layer.thickness)))
                self.targetTable.setCellWidget(row, 1+i*2, QLabel(str(layer)))
        self.targetTable.resizeColumnsToContents()
        self.targetTable.resizeRowsToContents()

    def update_reaction_table(self) -> None:
        self.reactionTable.setRowCount(len(self.spanc.reactions))
        self.reactionTable.setVerticalHeaderLabels(self.spanc.reactions.keys())
        for row, rxn in enumerate(self.spanc.reactions.values()):
            self.reactionTable.setItem(row, 0, QTableWidgetItem(str(rxn.targetMaterial)))
            self.reactionTable.setCellWidget(row, 1, QLabel(str(rxn)))
            self.reactionTable.setItem(row, 2, QTableWidgetItem(str(rxn.params.beamEnergy)))
            self.reactionTable.setItem(row, 3, QTableWidgetItem(str(rxn.params.magneticField)))
            self.reactionTable.setItem(row, 4, QTableWidgetItem(str(rxn.params.spsAngle / DEG2RAD)))
        self.reactionTable.resizeColumnsToContents()
        self.reactionTable.resizeRowsToContents()
        
    def update_calibration_table(self) -> None:
        self.calibrationTable.setRowCount(len(self.spanc.calibrations))
        self.calibrationTable.setVerticalHeaderLabels(self.spanc.calibrations.keys())
        for row, peak in enumerate(self.spanc.calibrations.values()):
            self.calibrationTable.setItem(row, 0, QTableWidgetItem(str(peak.peakID)))
            self.calibrationTable.setItem(row, 1, QTableWidgetItem(peak.rxnName))
            self.calibrationTable.setItem(row, 2, QTableWidgetItem(str(peak.position)))
            self.calibrationTable.setItem(row, 3, QTableWidgetItem(str(peak.positionErrStat)))
            self.calibrationTable.setItem(row, 4, QTableWidgetItem(str(peak.positionErrSys)))
            self.calibrationTable.setItem(row, 5, QTableWidgetItem(str(peak.rho)))
            self.calibrationTable.setItem(row, 6, QTableWidgetItem(str(peak.rhoErr)))
            self.calibrationTable.setItem(row, 7, QTableWidgetItem(str(peak.excitation)))
            self.calibrationTable.setItem(row, 8, QTableWidgetItem(str(peak.excitationErr)))
        self.calibrationTable.resizeColumnsToContents()
        self.calibrationTable.resizeRowsToContents()

    def update_output_table(self) -> None:
        self.outputTable.setRowCount(len(self.spanc.outputs))
        self.outputTable.setVerticalHeaderLabels(self.spanc.outputs.keys())
        for row, peak in enumerate(self.spanc.outputs.values()):
            self.outputTable.setItem(row, 0, QTableWidgetItem(str(peak.peakID)))
            self.outputTable.setItem(row, 1, QTableWidgetItem(peak.rxnName))
            self.outputTable.setItem(row, 2, QTableWidgetItem(str(peak.position)))
            self.outputTable.setItem(row, 3, QTableWidgetItem(str(peak.positionErrStat)))
            self.outputTable.setItem(row, 4, QTableWidgetItem(str(peak.positionErrSys)))
            self.outputTable.setItem(row, 5, QTableWidgetItem(str(peak.rho)))
            self.outputTable.setItem(row, 6, QTableWidgetItem(str(peak.rhoErr)))
            self.outputTable.setItem(row, 7, QTableWidgetItem(str(peak.excitation)))
            self.outputTable.setItem(row, 8, QTableWidgetItem(str(peak.excitationErr)))
            self.outputTable.setItem(row, 9, QTableWidgetItem(str(peak.positionFWHM)))
            self.outputTable.setItem(row, 10, QTableWidgetItem(str(peak.positionFWHMErr)))
            self.outputTable.setItem(row, 11, QTableWidgetItem(str(peak.excitationFWHM)))
            self.outputTable.setItem(row, 12, QTableWidgetItem(str(peak.excitationFWHMErr)))
        self.outputTable.resizeColumnsToContents()
        self.outputTable.resizeRowsToContents()

    def update_fit_order(self) -> None:
        self.fitOrderBox.setValue(self.spanc.fitter.polynomialOrder)

    #generate markdown text string and render in the text edit
    def update_fit_text(
        self,
        residuals: NDArray[np.float64],
        studentizedResiduals: NDArray[np.float64],
        energyResiduals: NDArray[np.float64],
        energyResidualErrors: NDArray[np.float64],
    ) -> None:
        markdownString = (
            f"# Fit Results\n\n"

            f"## Fit Statistics\n\n"
            f"- **Polynomial order:** {self.spanc.fitter.polynomialOrder}\n"
            f"- **Chi-square:** {self.spanc.fitter.get_chisquare():.3f}\n"
            f"- **NDF:** {self.spanc.fitter.get_ndf()}\n"
            f"- **Reduced chi-square:** "
            f"{self.spanc.fitter.get_reduced_chisquare():.3f}\n\n"

            f"## Calibration Parameters\n\n"
            f"**Parameter values ($a_0$ to $a_N$):**\n\n"
            f"```text\n"
            f"{np.array_str(self.spanc.fitter.get_parameters(), precision=3)}\n"
            f"```\n\n"

            f"**Parameter uncertainties ($u_{{a_0}}$ to $u_{{a_N}}$):**\n\n"
            f"```text\n"
            f"{np.array_str(self.spanc.fitter.get_parameter_errors(), precision=3)}\n"
            f"```\n\n"

            f"## Residual Diagnostics\n\n"
            f"**Rigidity residuals (cm, $x_0$ to $x_N$):**\n\n"
            f"```text\n"
            f"{np.array_str(residuals, precision=3)}\n"
            f"```\n\n"

            f"**Studentized residuals ($x_0$ to $x_N$):**\n\n"
            f"```text\n"
            f"{np.array_str(studentizedResiduals, precision=3)}\n"
            f"```\n\n"

            f"**Energy-equivalent residuals (keV, $x_0$ to $x_N$):**\n\n"
            f"```text\n"
            f"{np.array_str(energyResiduals, precision=2)}\n"
            f"```\n\n"

            f"**Energy-equivalent residual uncertainties "
            f"(keV, $x_0$ to $x_N$):**\n\n"
            f"```text\n"
            f"{np.array_str(energyResidualErrors, precision=2)}\n"
            f"```\n"
        )

        self.fitResultText.setMarkdown(markdownString)

def run_spanc_ui() :
    mpl.use("Qt5Agg")
    app = QApplication.instance()
    if not app:
        app = QApplication(sys.argv)
        app.setStyleSheet(load_stylesheet())
    window = SpancGUI()
    sys.exit(app.exec_())
