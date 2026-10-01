from .SPSReaction import *
from .SPSTarget import *
from .Fitter import *
from .data.NuclearData import *
from dataclasses import dataclass, field
import numpy as np
from enum import Enum

INVALID_PEAK_ID: int = -1
DEG2RAD: float = np.pi / 180.0

class PeakType(Enum):
    CALIBRATION = "Calibration"
    OUTPUT = "Output"

@dataclass
class Peak:
    excitation: float = 0.0 #MeV
    excitationErr: float = 0.0 #MeV
    position: float = 0.0 #arb
    positionErrStat: float = 0.0 #arb
    positionErrSys: float = 0.0 #arb
    rho: float = 0.0 #cm
    rhoErr: float = 0.0 #cm
    positionFWHM: float = 0.0
    positionFWHMErr: float = 0.0
    excitationFWHM: float = 0.0
    excitationFWHMErr: float = 0.0
    rxnName: str = ""
    peakID: int = INVALID_PEAK_ID

class Spanc:
    def __init__(self):
        self.targets: dict[str, SPSTarget] = {}
        self.reactions: dict[str, Reaction] = {}
        self.calibrations: dict[int, Peak] = {}
        self.outputs: dict[int, Peak] = {}
        self.fitter : Fitter = Fitter()
        self.energyFitters: dict[str, Fitter] = {}
        self.isFit: bool = False

    def ensure_energy_fitters(self) -> None:
        """Initialize attributes absent from SPANC files saved by older versions."""
        if not hasattr(self, "energyFitters"):
            self.energyFitters = {}
        for fitter in self.energyFitters.values():
            if not hasattr(fitter, "lastWarnings"):
                fitter.lastWarnings = []

    def set_fit_order(self, order: int) -> None:
        self.fitter.set_polynomial_order(order)

    #return fit data so that the data points can be drawn
    def fit(self) -> list[FitPoint]:
        fitData = [FitPoint(peak.position, peak.rho, np.sqrt(peak.positionErrStat**2.0 + peak.positionErrSys**2.0), peak.rhoErr) for peak in self.calibrations.values()]
        self.fitter.run(data=fitData)
        self.isFit = True
        return fitData

    def get_residuals(self) -> list[FitResidual]:
        return self.fitter.get_residuals()

    def has_current_rho_fit(self) -> bool:
        """Return whether the physical fit still matches all calibration data."""
        if (
            not self.isFit
            or self.fitter.fitResults is None
            or self.fitter.fitData is None
            or len(self.fitter.fitData) != len(self.calibrations)
        ):
            return False

        for point, peak in zip(
            self.fitter.fitData, self.calibrations.values()
        ):
            currentValues = (
                peak.position,
                peak.rho,
                np.hypot(peak.positionErrStat, peak.positionErrSys),
                peak.rhoErr,
            )
            fittedValues = (point.x, point.y, point.xError, point.yError)
            if not np.allclose(
                fittedValues,
                currentValues,
                rtol=1.0e-12,
                atol=1.0e-12,
                equal_nan=True,
            ):
                return False
        return True

    def get_calibrations_for_reaction(self, rxnName: str) -> list[Peak]:
        return [
            peak
            for peak in self.calibrations.values()
            if peak.rxnName == rxnName
        ]

    def get_energy_fitter(self, rxnName: str) -> Fitter:
        self.ensure_energy_fitters()
        if rxnName not in self.reactions:
            raise ValueError(f"Unknown reaction '{rxnName}'.")
        if rxnName not in self.energyFitters:
            self.energyFitters[rxnName] = Fitter(order=1)
        return self.energyFitters[rxnName]

    def set_energy_fit_order(self, rxnName: str, order: int) -> None:
        fitter = self.get_energy_fitter(rxnName)
        if fitter.polynomialOrder != order:
            fitter.invalidate()
            fitter.set_polynomial_order(order)

    def invalidate_energy_fit(self, rxnName: str) -> None:
        self.ensure_energy_fitters()
        if rxnName in self.energyFitters:
            self.energyFitters[rxnName].invalidate()

    def has_energy_fit(self, rxnName: str) -> bool:
        self.ensure_energy_fitters()
        return (
            rxnName in self.energyFitters
            and self.energyFitters[rxnName].fitResults is not None
            and self.energyFitters[rxnName].function is not None
        )

    def fit_energy(self, rxnName: str) -> list[FitPoint]:
        """Fit a reaction-specific focal-plane position to energy calibration."""
        fitter = self.get_energy_fitter(rxnName)
        peaks = self.get_calibrations_for_reaction(rxnName)
        minimumPoints = fitter.polynomialOrder + 2
        if len(peaks) < minimumPoints:
            raise ValueError(
                f"Reaction {rxnName} has {len(peaks)} calibration point(s), but "
                f"an order-{fitter.polynomialOrder} fit requires at least "
                f"{minimumPoints}."
            )

        fitData = [
            FitPoint(
                peak.position,
                peak.excitation,
                np.sqrt(
                    peak.positionErrStat**2.0
                    + peak.positionErrSys**2.0
                ),
                peak.excitationErr,
            )
            for peak in peaks
        ]
        fitter.run(data=fitData, allowMissingErrors=True)
        return fitData

    @staticmethod
    def format_polynomial(
        parameters: np.ndarray,
        scale: float = 1.0,
        variable: str = "x",
    ) -> str:
        """Format ascending polynomial coefficients as a copyable equation."""
        terms = []
        for order, rawCoefficient in enumerate(parameters):
            coefficient = float(rawCoefficient) * scale
            if coefficient == 0.0 and len(parameters) > 1:
                continue
            magnitude = f"{abs(coefficient):.10g}"
            if order == 0:
                body = magnitude
            elif order == 1:
                body = f"{magnitude} * {variable}"
            else:
                body = f"{magnitude} * {variable}^{order}"

            if not terms:
                terms.append(f"-{body}" if coefficient < 0.0 else body)
            else:
                terms.append(
                    f" {'-' if coefficient < 0.0 else '+'} {body}"
                )
        return "".join(terms) if terms else "0"

    def get_energy_equations(self, rxnName: str) -> tuple[str, str]:
        if not self.has_energy_fit(rxnName):
            raise ValueError(f"Reaction {rxnName} does not have an energy fit.")
        parameters = self.get_energy_fitter(rxnName).get_parameters()
        equationMeV = f"Ex(MeV) = {self.format_polynomial(parameters)}"
        equationKeV = f"Ex(keV) = {self.format_polynomial(parameters, 1000.0)}"
        return equationMeV, equationKeV

    def evaluate_energy(self, rxnName: str, position: float) -> float:
        if not self.has_energy_fit(rxnName):
            raise ValueError(f"Reaction {rxnName} does not have an energy fit.")
        return float(self.get_energy_fitter(rxnName).evaluate(position))

    def get_energy_sensitivity(self, rxnName: str, position: float) -> float:
        """Return dEx/dx at position in keV/mm."""
        if not self.has_energy_fit(rxnName):
            raise ValueError(f"Reaction {rxnName} does not have an energy fit.")
        return float(
            self.get_energy_fitter(rxnName).evaluate_derivative(position)
            * 1000.0
        )

    def get_direct_energy_residuals(self, rxnName: str) -> np.ndarray:
        """Return assigned-minus-fitted direct energy residuals in keV."""
        if not self.has_energy_fit(rxnName):
            return np.array([])
        fitter = self.get_energy_fitter(rxnName)
        return np.asarray(
            [residual.residual * 1000.0 for residual in fitter.get_residuals()]
        )

    def get_energy_residuals(self) -> tuple[np.ndarray, np.ndarray]:
        """Convert calibration residuals from rho space to energy space.

        Returns
        -------
        energyResiduals
            Assigned excitation energy minus the excitation energy predicted
            by the fitted calibration, in keV.
        energyResidualErrors
            Rho-space residual uncertainties converted to equivalent
            excitation-energy uncertainties, in keV.

        Notes
        -----
        This conversion is performed only after the ODR fit, so it does not
        alter the fit or its statistical weighting.
        """
        if not self.has_current_rho_fit():
            return np.array([]), np.array([])

        rhoResiduals = self.get_residuals()
        energyResiduals = np.empty(len(rhoResiduals))
        energyResidualErrors = np.empty(len(rhoResiduals))

        for i, (peak, residual) in enumerate(
            zip(self.calibrations.values(), rhoResiduals)
        ):
            rxn = self.reactions[peak.rxnName]

            # Rho predicted by the fitted calibration at this position.
            rhoFit = self.fitter.evaluate(peak.position)
            exFit = rxn.calculate_excitation(rhoFit)

            # Assigned energy minus fitted energy, reported in keV.
            energyResiduals[i] = (peak.excitation - exFit) * 1000.0

            # Convert the combined rho-residual uncertainty through the full
            # reaction kinematics using a symmetric interval about rhoFit.
            exPlus = rxn.calculate_excitation(
                rhoFit + residual.residualError
            )
            exMinus = rxn.calculate_excitation(
                rhoFit - residual.residualError
            )
            energyResidualErrors[i] = (
                0.5 * abs(exPlus - exMinus) * 1000.0
            )

        return energyResiduals, energyResidualErrors

    def add_target(self, targName: str, layers: list[TargetLayer]) -> None:
        self.targets[targName] = SPSTarget(layers, name=targName)

    def add_reaction(self, params: RxnParameters, targetName: str) -> None:
        if targetName not in self.targets:
            print("Cannot create reaction with non-existant target ", targetName)
            return
        key = f"Rxn{len(self.reactions)}"
        params.spsAngle *= DEG2RAD
        rxn = Reaction(params, target=self.targets[targetName])
        self.reactions[key] = rxn

    def update_reaction_parameters(self, beamEnergy: float, spsAngle: float, magneticField: float, rxnName: str):
        if rxnName in self.reactions:
            rxn = self.reactions[rxnName]
            rxn.params.beamEnergy = beamEnergy
            rxn.params.spsAngle = spsAngle * DEG2RAD
            rxn.params.magneticField = magneticField

    def add_calibration(self, data: Peak) -> None:
        if data.rxnName in self.reactions:
            oldRxnName = None
            if data.peakID in self.calibrations:
                oldRxnName = self.calibrations[data.peakID].rxnName
            rxn = self.reactions[data.rxnName]
            data.rho = rxn.convert_ejectile_KE_2_rho(rxn.calculate_ejectile_KE(data.excitation))
            data.rhoErr = np.abs(rxn.convert_ejectile_KE_2_rho(rxn.calculate_ejectile_KE(data.excitation + data.excitationErr)) - data.rho)
            if data.peakID == INVALID_PEAK_ID:
                data.peakID = len(self.calibrations)
            self.calibrations[data.peakID] = data
            if oldRxnName is not None:
                self.invalidate_energy_fit(oldRxnName)
            self.invalidate_energy_fit(data.rxnName)
        return
    
    def remove_calibration(self, data: Peak) -> bool:
        if data.peakID not in self.calibrations.keys():
            return False
        rxnName = self.calibrations[data.peakID].rxnName
        self.calibrations.pop(data.peakID)
        self.invalidate_energy_fit(rxnName)
        return True

    def add_output(self, data: Peak) -> None:
        if data.peakID == INVALID_PEAK_ID:
            data.peakID = len(self.outputs)
        self.outputs[data.peakID] = data
        return

    def calculate_output_urho(self, peak: Peak) -> float:
        urho = 0.0
        paramErrors = self.fitter.get_parameter_errors()
        for i, paramErr in enumerate(paramErrors):
            urho += (self.fitter.evaluate_param_derivative(peak.position, i) * paramErr)**2.0
        urho += (self.fitter.evaluate_derivative(peak.position)*np.sqrt(peak.positionErrStat**2.0 + peak.positionErrSys**2.0))**2.0
        return np.sqrt(urho)

    def calculate_outputs(self) -> None:
        if self.isFit == False:
            return

        for output in self.outputs.values():
            rxn = self.reactions[output.rxnName]
            output.rho = self.fitter.evaluate(output.position)
            output.rhoErr = self.calculate_output_urho(output)
            output.excitation = rxn.calculate_excitation(output.rho)
            output.excitationErr = np.abs(rxn.calculate_excitation(output.rho + output.rhoErr) - output.excitation)
			
            if output.positionFWHM == 0:
                output.excitationFWHM = 0
                output.excitationFWHMErr = 0
            else:
                rhoLo = self.fitter.evaluate(output.position - output.positionFWHM * 0.5)
                rhoHi = self.fitter.evaluate(output.position + output.positionFWHM * 0.5)
                exLo = rxn.calculate_excitation(rhoLo)
                exHi = rxn.calculate_excitation(rhoHi)
                output.excitationFWHM = abs(exHi - exLo)
                output.excitationFWHMErr = output.positionFWHMErr/output.positionFWHM*output.excitationFWHM

    def calculate_calibrations(self) -> None:
        for calibration in self.calibrations.values():
            rxn = self.reactions[calibration.rxnName]
            calibration.rho = rxn.convert_ejectile_KE_2_rho(rxn.calculate_ejectile_KE(calibration.excitation))
            calibration.rhoErr = np.abs(rxn.convert_ejectile_KE_2_rho(rxn.calculate_ejectile_KE(calibration.excitation + calibration.excitationErr)) - calibration.rho)
