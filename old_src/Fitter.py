import numpy as np
from numpy.typing import NDArray
from numpy.polynomial import Polynomial
from scipy import odr
from dataclasses import dataclass
from typing import Optional

INVALID_FIT_RESULT = np.inf
INVALID_NDF = -1


@dataclass
class FitPoint:
    x: float = 0.0
    y: float = 0.0
    xError: float = 0.0
    yError: float = 0.0


def convert_fit_points_to_arrays(
    data: list[FitPoint],
) -> tuple[
    NDArray[np.float64], NDArray[np.float64], NDArray[np.float64], NDArray[np.float64]
]:
    xArray = np.empty(len(data))
    yArray = np.empty(len(data))
    xErrorArray = np.empty(len(data))
    yErrorArray = np.empty(len(data))
    for index, point in enumerate(data):
        xArray[index] = point.x
        yArray[index] = point.y
        xErrorArray[index] = point.xError
        yErrorArray[index] = point.yError
    return xArray, yArray, xErrorArray, yErrorArray


@dataclass
class FitResidual:
    x: float = 0.0
    residual: float = 0.0
    residualError: float = 0.0
    studentizedResidual: float = 0.0


def convert_resid_points_to_arrays(
    data: list[FitResidual],
) -> tuple[
    NDArray[np.float64], NDArray[np.float64], NDArray[np.float64], NDArray[np.float64]
]:
    xArray = np.empty(len(data))
    residArray = np.empty(len(data))
    residErrorArray = np.empty(len(data))
    studentResidArray = np.empty(len(data))
    for index, point in enumerate(data):
        xArray[index] = point.x
        residArray[index] = point.residual
        residErrorArray[index] = point.residualError
        studentResidArray[index] = point.studentizedResidual
    return xArray, residArray, residErrorArray, studentResidArray


class Fitter:
    def __init__(self, order: int = 1):
        self.polynomialOrder: int = order
        self.fitResults: Optional[odr.Output] = None
        self.fitData: Optional[list[FitPoint]] = None
        self.function: Optional[Polynomial] = None

    def set_polynomial_order(self, order: int) -> None:
        self.polynomialOrder = order

    def run(self, data: list[FitPoint] = None) -> None:
        if data is not None:
            self.fitData = data

        if self.fitData is not None:
            xArray, yArray, xErrorArray, yErrorArray = convert_fit_points_to_arrays(
                self.fitData
            )
            modelData = odr.RealData(xArray, y=yArray, sx=xErrorArray, sy=yErrorArray)
            model = odr.polynomial(self.polynomialOrder)
            self.fitResults = odr.ODR(modelData, model).run()
            self.function = Polynomial(self.fitResults.beta)
        else:
            print("Cannot run fitter without setting data to be fit!")

    def get_parameters(self) -> NDArray[np.float64]:
        if self.fitResults is not None:
            return self.fitResults.beta
        return np.array({INVALID_FIT_RESULT})

    def get_parameter_errors(self) -> NDArray[np.float64]:
        if self.fitResults is not None:
            return self.fitResults.sd_beta
        return np.array({INVALID_FIT_RESULT})

    def get_ndf(self) -> int:
        """Return the number of degrees of freedom for the polynomial fit.

        A polynomial of order ``m`` has ``m + 1`` fitted coefficients, so
        NDF = N_data - (m + 1).
        """
        if self.fitResults is not None and self.fitData is not None:
            n_parameters = len(self.fitResults.beta)
            return len(self.fitData) - n_parameters
        return INVALID_NDF

    def evaluate(self, x: float) -> float:
        if self.function is not None:
            return self.function(x)
        return INVALID_FIT_RESULT

    def evaluate_derivative(self, x: float) -> float:
        if self.function is not None:
            return self.function.deriv()(x)
        return INVALID_FIT_RESULT

    def evaluate_param_derivative(self, x: float, index: int) -> float:
        if self.fitResults is not None and len(self.fitResults.beta) > index:
            return x**index
        return INVALID_FIT_RESULT

    def get_chisquare(self) -> float:
        """Return ODR's weighted sum of squares (chi-square-like statistic)."""
        if self.fitResults is None:
            return INVALID_FIT_RESULT
        return float(self.fitResults.sum_square)

    def get_reduced_chisquare(self) -> float:
        """Return the weighted sum of squares per degree of freedom."""
        if self.fitResults is None:
            return INVALID_FIT_RESULT

        ndf = self.get_ndf()
        if ndf <= 0:
            return INVALID_FIT_RESULT

        return self.get_chisquare() / ndf

    def get_residuals(self) -> list[FitResidual]:
        if self.fitData is None or self.fitResults is None:
            return []

        # For r = y - f(x), propagate the measurement uncertainties in both
        # coordinates into the vertical residual direction:
        #
        #   sigma_r^2 = sigma_y^2 + [f'(x) sigma_x]^2
        #
        # This matches the x/y uncertainty information supplied to the ODR fit.
        fitResiduals = []
        for point in self.fitData:
            residual = point.y - self.evaluate(point.x)
            residualError = np.sqrt(
                point.yError**2.0
                + (self.evaluate_derivative(point.x) * point.xError) ** 2.0
            )
            fitResiduals.append(
                FitResidual(point.x, residual, residualError, 0.0)
            )

        ndf = self.get_ndf()
        if ndf <= 0:
            return fitResiduals

        # Studentized residuals are used here only as a diagnostic for the
        # polynomial calibration.  Compute the leverage from the polynomial
        # design matrix, H = X (X^T X)^(-1) X^T, and use the usual internally
        # studentized residual definition r_i / [s sqrt(1 - h_ii)].
        #
        # Note: the fitted coefficients themselves come from ODR (which can
        # include uncertainties in both x and y), so these studentized
        # residuals should be interpreted as an approximate diagnostic rather
        # than as an exact ODR goodness-of-fit statistic.
        x = np.asarray([point.x for point in self.fitData], dtype=float)
        residuals = np.asarray([resid.residual for resid in fitResiduals], dtype=float)

        design = np.vander(x, N=self.polynomialOrder + 1, increasing=True)
        hat = design @ np.linalg.pinv(design.T @ design) @ design.T
        leverage = np.clip(np.diag(hat), 0.0, 1.0)

        rss = float(np.sum(residuals**2))
        residual_std = np.sqrt(rss / ndf)

        if residual_std == 0.0:
            return fitResiduals

        for resid, h_ii in zip(fitResiduals, leverage):
            denom = residual_std * np.sqrt(max(1.0 - h_ii, np.finfo(float).eps))
            resid.studentizedResidual = resid.residual / denom

        return fitResiduals
