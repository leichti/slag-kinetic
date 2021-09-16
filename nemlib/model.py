import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import least_squares, fsolve
import statistics

class KineticModel:
    def __init__(self, trial, dimensions, t0_0=True):
        """
        Base class for numerical modelling of ODEs, based on solve_ivp().
        Child Classes must only implement the properties
        KineticModel.initial_guess, KineticModel.bounds, KineticModel.used_model_compounds
        and the system of differential equations representing the system's reactions
        in KineticModel.model
        It provides the following:
            - A wrapper to fit the model to a given SlagReductionTrial.
            - A user-friendly methods to get mole or mole concentration over time for
              given model parameters or optimized model parameters.
            - A method to get the deviation between experimental and model values
              vs the experimental values to easily generate deviation plots
            - A method to get the absolute reaction rate for given
              experimental values. The way it works is to calculate y-values
              over time based on KineticModel.optimized_parameters and subsequent
              building the first derivative. Then it returns both lists, the y-values and first derivative
              of the y-values. For many equations it is possible to use the
              differential representation in the KineticModel.model() method, however,
              in some cases, the reaction rate for a given compound depends on multiple
              compounds, rule out this easier approach (e.g. reaction rate of FeO
              in the ZnOFeModel is a function of x_FeO and x_ZnO)
            - A method KineticModel.shift_to_y0(n) to shift the time axis to get y(t=0)=n.
              This is helpful to compare multiple experiments in a single  plot, all having one same
              y0 value at t0.

        Parameters
        ----------
        trial : SlagReductionTrial
            SlagReductionTrial representing the (absolute) mol course over time
            of the experiment.
        dimensions : TrialDimensions
            Representing the dimensions of the experiment. Directly
            including those allows to directly fit mass transfer coefficients
            in the KineticModel.model() method.
        t0_0 : Option to shift the experimental data to start with the first
            sample at t0=0.

        """
        self.trial = trial
        self.optimized_parameters = []
        self.dimension = dimensions

        if t0_0 is True:
            self.trial.t_shift(-self.trial.t_min)

    @property
    def initial_guess(self):
        """

        Returns
        -------
        initial_guess : list or tuple or array-like
                Array-like list of values with a length that equals
                the number of independent variables plus number of
                arguments (e.g. ya, yb, yc.., k1, k2, .. kn)
        """
        raise NotImplementedError("Please implement initial_guess property")

    @property
    def used_model_compounds(self):
        """

        Returns
        -------
        number : int
            number of model compounds used in the reaction system
        """
        raise NotImplementedError("Please specify how many model compounds are used"
                                  "in your implementation")

    @property
    def bounds(self):
        """
        Returns
        -------
        bounds : tuple
            bounds for the parameters in the least square's "trf" algorithm
            If the child class doesn't include a method, default values are
            -inf to +inf for all parameters.

        """
        lower_bounds = []
        upper_bounds = []

        for _ in self.initial_guess:
            lower_bounds.append(-np.inf)
            upper_bounds.append(np.inf)

        return tuple([lower_bounds, upper_bounds])

    def model(self, t, y0, *args):
        """
        The abstract model method representing the differential kinetic model.
        Current implementation doesn't allow vectorized calculations

        Parameters
        ----------
        t : float
            time at which the reaction rates are calculated
        y0 : float or list
            value of the independent variables at time t. If the system
            has only one independent variable, the type of y is float,
            otherwise list
        args : list
            the kinetic parameters for the model

        Returns
        -------
        rates : float or array-like
            returns the reaction rate (dy/dt) or reaction rates
            ([dy1/dt, dy2/dt, ..., dyi/dt) for the given parameters
        """
        raise NotImplementedError("Method model not implemented")

    @property
    def fit_report(self):
        """

        Returns
        -------
        square_sum : float
            sum of squared deviation between experimental and model values
        """
        dev = self.fit_function(self.optimized_parameters)
        return f"{statistics.mean(dev):.2f} ± {statistics.stdev(dev):.2f}"

    def rr_v_y(self, which=0, t=None):
        """

        Parameters
        ----------
        t : None or list
            times at which the reaction rate and corresponding
            y-value should be determined
        which : int
            The integer key of the model compound.
            e.g. the first compound has the key 0, the second the key 1 and so on

        Returns
        -------
        y_pct, rr : np.array, np.array
            number pairs of relative y-values and the corresponding reaction
            rate at this y-values.

        """
        if t is None:
            t = np.arange(*self.trial.t_min_max, 0.1)
        y_pct = self.y_pct(t, *self.optimized_parameters)[which]
        rr = np.diff(y_pct) / np.diff(t)
        return y_pct[0:-1], -rr

    def deviation(self, which, pct=True):
        """
        Calculates and returns the absolute deviation (mol) between
        the model with KineticModel.optimized_parameters and the corresponding
        experimental values from the SlagReductionTrial

        Parameters
        ----------
        which : int
            The model compound (y_which) for which the reaction rate
            of the should be calculated. STARTS WITH 0!

        pct : bool
            if false: returns absolute deviation in mol
            if true: returns deviation in mol-%

        Returns
        -------
        times : list
            float time representatives at which samples where taken.
            e.g. seconds, minutes, hours from start

        deviation : np.array
            the deviation between model and experimental values

        """
        values = self.trial.model_compounds.x(which, pct=pct)
        times = self.trial.t

        if pct is True:
            model_values = self.y_pct(times, *self.optimized_parameters)[which]
        else:
            model_values = self.y(times, *self.optimized_parameters)[which]

        return times, np.subtract(model_values, values)

    def y(self, t, *args):
        """
        Based on solve_ivp. Solves the system of differential
        equations from the model() method

        Parameters
        ----------
        t : list
            time points at which y values should be returned
        args: list of arguments
            contains initial_conditions followed by the model parameters

        Returns
        -------
        y : np.array
            y values in mol
        """
        ivp_result = solve_ivp(self.model, self.t_span(t),
                               dense_output=True,
                               y0=self.y0(args),
                               args=self.args(args),
                               method="RK45")
        model_values = ivp_result.sol(t)

        # fix that moles can for some models be negative mathematically.
        for i, value in enumerate(model_values):
            model_values[i] = np.where(value < 0, 0, value)

        return model_values

    def y_pct(self, *args):
        """
        Wrapper around the y-method to get mol-pct instead of absolute mol

        Parameters
        ----------
        args : arguments

        Returns
        -------
        y_pct : np.array
            y values in mol-%
        """
        mol = self.y(*args)
        for i, v in enumerate(mol):
            mol[i] = np.divide(v, (self.trial.inert_moles + v)) * 100

        return mol

    def fit(self):
        """
        Fits model y() to the experimental values. Based on least_squares.
        In first step a rough estimation using the trf region algorithm,
        The rough estimated parameters from the first step
        are used as starting point in the second fit using the lm algorithm.

        Returns
        -------
        parameters : list
            list with the optimized parameters

        """
        fit_result = least_squares(self.fit_function,
                                   x0=self.initial_guess,
                                   bounds=self.bounds,
                                   method="trf")
        fit_result = least_squares(self.fit_function,
                                   x0=fit_result["x"],
                                   method="lm")
        self.optimized_parameters = fit_result["x"]
        return fit_result["x"]

    def fit_function(self, parameters):
        """
        Calculates a list holding the deviation between model and
        experimental values.

        Parameters
        ----------
        parameters : list
            list of model parameters

        Returns
        -------
        deviation : np.array
            list of deviation values (model - experimental)

        """

        model_values = self.y_pct(self.trial[self.trial.time_col], *parameters)
        dev2 = np.array([])

        for i, model_value in enumerate(model_values):
            compound = self.trial.model_compounds.name_x(i)
            dev = np.subtract(model_value,
                              self.trial.y(compound, pct=True))
            dev2 = np.append(dev2, dev)

        return dev2

    def y0(self, args=None):
        """
        Helper function. Extracts the y0_values from the parameters
        (e.g. y0, y1, .., yn, k1, k2, ..., ki extracts y0-yn)

        Parameters
        ----------
        args : list
            arguments for the model. Usually the KineticModel.optimized_parameters

        Returns
        -------
        y0_values : list or float
            y0 values from the arguments list

        """
        if args is None:
            args = self.optimized_parameters

        return args[0: self.used_model_compounds]

    def args(self, args=None):
        """
        Helper function. Extracts the y0_values from the parameters
        (e.g. y0, y1, .., yn, k1, k2, ..., ki extracts k1-ki)

        Parameters
        ----------
        args : list
            arguments for the model. Usually the KineticModel.optimized_parameters

        Returns
        -------
        y0_values : list or float
            model parameter values from the arguments list

       """
        if args is None:
            args = self.optimized_parameters

        return args[self.used_model_compounds:]

    @staticmethod
    def t_span(t):
        """
        Helper function for the y method() to determine the time frame in which
        the numerical solution should be solved.

        Parameters
        ----------
        t : list
            a list with time points

        Returns
        -------
        min : float
            the minimal value. Is always 0 or less.
        max: float
            the maximum value

        """
        if min_t := min(t) > 0:
            min_t = 0

        return min_t, max(t)

    def shift_to_y0(self, y, which=0):
        """
        Shifts the model and experimental time values to get f(t=0)=y.
        This is helpful if multiple experiments are plotted in a single plot.
        Shifting requires a fitted model. Then fsolve is used to find
        t where f(t) = y. The found time value is subtracted from the
        experimental time values. Finally, the model is fitted again to the
        time-shifted experimental values. This approach is required for
        models with more than 1 model compound
        (e.g. ZnO and FeO, ZnO_0=y is given, but FeO_0 must be fitted
         again after shifting)


        Parameters
        ----------
        y : float
            wanted y0 value
        which : int
            the model compound on which the shifting is based on.

        Returns
        -------
        result : list
            list of time values
        """
        if len(self.optimized_parameters) == 0:
            self.fit()

        result = fsolve(self.time_shift_function, np.array([0]), args=(y, which))

        if len(result) == 0:
            raise ValueError("No solution found.")

        self.trial.t_shift(-result[0])
        self.fit()

        return result

    def time_shift_function(self, t, y0_wanted, which):
        """
        Helper function for scipy.optimize.fsolve in shift_to_y0
        Calculates the y_pct value for a given t and subtracts
        the wanted y value.

        Parameters
        ----------
        t : float
        y0_wanted : float
        which : int
            the model compound on which the shifting is based on.
            (e.g. 0 could be ZnO, 1 could be FeO).
            BE AWARE, multiple solutions are possible
            (e.g. the FeO in the ZnFeOModel)

        Returns
        -------
        deviation : float
            deviation between model value and wanted value

        """
        dev = self.y_pct(t, *self.optimized_parameters)[which] - y0_wanted
        return dev


class FirstOrderModel(KineticModel):

    def __init__(self, *args):
        """
        A simple first order model for one model compound
        dy/dt = -k*y


        """
        super().__init__(*args)
        self.k = 0

    def __str__(self):
        return r"$\dfrac{dy_place}{dt}=-k_{reduction} \cdot A \cdot x_{y_place}$"

    def label(self, which=0, simple=True):
        oxide = self.trial.model_compounds.name_x(which)

        if simple is True:
            return r"$\dfrac{d %s}{dt}= -k_2 \cdot x_{%s}$" % (oxide, oxide)

        area = self.dimension.area
        k = self.optimized_parameters[1:][0]
        return r"$\dfrac{d%s}{dt}=-%.2f \cdot %.0f \cdot x_{%s}$" \
               % (oxide, k, area, oxide)

    @property
    def used_model_compounds(self):
        return 1

    @property
    def bounds(self):
        return tuple([(0, 0), (np.inf, np.inf)])

    @property
    def initial_guess(self):
        return tuple([self.trial.model_compounds.first[0], 0])

    def model(self, t, y, *parameters):
        k_reduction = parameters[0]  # type: float
        total_moles = y + self.trial.inert_moles
        dy_dt = -k_reduction * self.dimension.area * np.divide(y, total_moles) / 1000
        return dy_dt


class ZnOCCrucibleModel(KineticModel):
    def __init__(self, *args, n=1):
        super().__init__(*args)
        self.n = n

    def label(self, which=0, simple=True):
        oxide = self.trial.model_compounds.name_x(which)

        if simple is True:
            return r"$\dfrac{d %s}{dt}= -k_2 \cdot x_{%s} " \
                   r"- k_3 $" % (oxide, oxide)

        area = self.dimension.area
        k_zno_c, k_oxidation = self.optimized_parameters[1:]
        if k_oxidation < 0:
            sign = "-"
            k_oxidation = -k_oxidation
        else:
            sign = "+"
        return r"$\dfrac{d%s}{dt}=-%.2f \cdot %.0f \cdot x_{%s} %s %.2f$" \
               % (oxide, k_zno_c, area, oxide, sign, k_oxidation)

    @property
    def bounds(self):
        return (0, 0, -np.inf), (np.inf, np.inf, np.inf)

    @property
    def initial_guess(self):
        return tuple([self.trial.model_compounds.first[0], 0.01, -0.1])

    @property
    def used_model_compounds(self):
        return 1

    def model(self, t, y, *parameters):
        k_zno_c, k_oxidation = parameters
        zno_moles = y

        area = self.dimension.area
        total_moles = self.trial.inert_moles + zno_moles

        d_moles_dt = -k_zno_c * area * np.power(np.divide(zno_moles, total_moles), self.n) / 1000 + \
                     k_oxidation / 1000

        return d_moles_dt


class NickelModel(KineticModel):
    def __init__(self, *args, n=1):
        super().__init__(*args)
        self.n = n

    def label(self, which=0, simple=True):
        oxide = self.trial.model_compounds.name_x(which)

        if simple is True:
            return r"$\dfrac{d %s}{dt}= -k_1 \cdot x_{%s}^2 - k_2 \cdot x_{%s} - k_3$" \
                   % (oxide, oxide,oxide)

        area_metal_c = self.dimension.bottom
        area_crucible_c = self.dimension.shell
        k2_reduction, k_reduction, k_oxidation, *_ = self.optimized_parameters[1:]

        return r"$\dfrac{d%s}{dt}\;=" \
               r" -%.2f \cdot %.0f \cdot x_{%s}^2 $" \
                " \n" \
               r"$ + %.2f \cdot %.0f \cdot x_{%s}  - %.2f$" \
               % (oxide, k2_reduction, area_metal_c, oxide, k_reduction,
                  area_crucible_c, oxide, k_oxidation)

    @property
    def bounds(self):
        return (0, 0, -np.inf, -np.inf), (np.inf, np.inf, np.inf, np.inf)

    @property
    def initial_guess(self):
        return tuple([self.trial.model_compounds.first[0], 0, 0, 0])

    @property
    def used_model_compounds(self):
        return 1

    def model(self, t, y, *parameters):
        k2_zno_c, k_zno_c, k_oxidation = parameters

        area_metal_c = self.dimension.bottom
        area_crucible_c = self.dimension.shell
        total_moles = self.trial.inert_moles + y
        x1 = np.divide(y, total_moles)

        d_moles_dt = - k2_zno_c * area_metal_c * np.power(x1, 2) / 1000 \
                     - k_zno_c * area_crucible_c * x1 / 1000 \
                     + k_oxidation / 1000

        # must return a single value
        return d_moles_dt


class ZnOFeModel(KineticModel):

    def __init__(self, *args):
        super().__init__(*args)

    @property
    def initial_guess(self):
        return (self.trial.model_compounds.first[0],
                self.trial.model_compounds.second[0],
                0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1)

    @property
    def bounds(self):
        return ((0, 0, 0, 0, 0, 0, 0, 0, 0),
                (np.inf, np.inf, np.inf, np.inf, np.inf, np.inf, np.inf, np.inf, np.inf))

    @property
    def used_model_compounds(self):
        return 2

    def label(self, which=0, simple=True):
        oxide = self.trial.model_compounds.name_x(which)

        if simple is True:
            return r"$\dfrac{d %s}{dt}=$ ZnOFeModel" % oxide

        t, y, k_zno_c_metal, k_zno_c, k_zno_fe, k_zn_oxidation, \
            k_feo_c_metal, k_feo_c, k_fe_oxidation = self.optimized_parameters

        shell_area = self.dimension.shell
        bottom_area = self.dimension.bottom
        if which == 0:
            return r"$\dfrac{d %s}{dt}=$ " \
                   r"$-(%.2f+%.2f) \cdot %.0f \cdot x_{%s}^2 $" \
                   r"$-%.2f \cdot %.0f \cdot x_{%s} $" \
                   r"$-%.2f$" \
                   % (oxide,
                      k_zno_c_metal, k_zno_fe, bottom_area, oxide,
                      k_zno_c, shell_area, oxide,
                      k_zn_oxidation)
        else:

            oxide2 = self.trial.model_compounds.name_first
            return r"$\dfrac{d %s}{dt}= " \
                   r"-(%.2f) \cdot %.0f \cdot x_{%s}^2" \
                   r"-%.2f \cdot %.0f \cdot x_{%s} " \
                   r"-%.2f \cdot %.0f \cdot x_{%s}^2 " \
                   r"-%.2f$" \
                   % (oxide,
                      k_feo_c_metal, bottom_area, oxide,
                      k_feo_c, shell_area, oxide,
                      k_zno_fe, bottom_area, oxide2,
                      k_fe_oxidation)

    def model(self, t, y, *parameters):
        k_zno_c_metal, k_zno_c, k_zno_fe, k_zn_oxidation, \
            k_feo_c_metal, k_feo_c, k_fe_oxidation = parameters
        zno_moles, feo_moles = y
        bottom = self.dimension.bottom
        shell = self.dimension.shell
        total_moles = self.trial.inert_moles + zno_moles
        x_zno = np.divide(zno_moles, total_moles)
        x_feo = np.divide(feo_moles, total_moles)
        d_zno_dt = - (k_zno_c_metal + k_zno_fe) * bottom * np.power(x_zno, 2) / 1000 \
                   - k_zno_c * shell * x_zno / 1000 \
                   - k_zn_oxidation / 1000
        # must return a single value
        d_feo_dt = - k_feo_c_metal * bottom * np.power(x_feo, 2) / 10000 \
                   - k_feo_c * shell * x_feo / 1000 \
                   + k_zno_fe * bottom * np.power(x_zno, 2) / 1000 \
                   + k_fe_oxidation / 1000

        return d_zno_dt, d_feo_dt
