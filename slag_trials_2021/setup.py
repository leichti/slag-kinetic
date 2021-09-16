import copy

import numpy as np
from nemlib.analysis import PhaseConstructor, ElementalAnalysesFile, InfoOrganizer
from nemlib.model import ZnOCCrucibleModel, ZnOFeModel, KineticModel
from nemlib.modelhelper import TrialDimensions, TrialAssembler, \
    ModelCompounds
from nemlib.trial import SlagReductionTrial, VxPySelector
import matplotlib.pyplot as plt


class MetalBathSeries:

    def __init__(self):
        self.slag_density = 2500
        ##############################
        # Load and organize analyses #
        ##############################
        phase_constructor = PhaseConstructor(phases=["CaF2", "CaO", "SiO2",
                                                     "Al2O3", "MgO", "Na2O",
                                                     "K2O", "FeO", "ZnO"],
                                             ignore=["O", "C"])
        self.phased_data = ElementalAnalysesFile("data/sem.xlsx").phased(phase_constructor)
        time_data = InfoOrganizer("data/time_table.xlsx", "sample", ["time"])
        self.phased_data.add_info(time_data)
        trial_info = InfoOrganizer("data/trial_info.xlsx", "trial", ["inert mass", "zno mass", "metal mass", "setup"])

        inert_compounds = ["CaO", "SiO2", "Al2O3", "MgO", "CaF2", "K2O", "Na2O"]

        rate_constants = []

        model_compounds = ModelCompounds(["ZnO", "FeO"])
        self.assembler = TrialAssembler(trial_info=trial_info, inert_compounds=inert_compounds,
                                        dimensions=TrialDimensions, slag_density=2.5,
                                        model=ZnOCCrucibleModel, model_compounds=model_compounds)


class Plot:
    def __init__(self, setup):
        self.phased_data = setup.phased_data
        self.assembler = setup.assembler
        self.setup = setup

class SummaryPlot(Plot):
    def __init__(self, setup, titles=[]):
        super(SummaryPlot, self).__init__(setup)

        self.fig, self.ax = plt.subplots(ncols=2)
        self.titles = titles
        self.trials = []

    def add(self, name, droplets=[], model=None, add_scatter=False, which=0):
        trial = SlagReductionTrial(VxPySelector(name)(self.phased_data))

        if model is not None:
            self.assembler.model = model

        trial = self.assembler.do(trial, exclude=droplets)
        trial.model.shift_to_y0(12)

        self.ax[0].plot(*trial.model_y(pct=True, t=np.arange(0, 200, 0.1)), label=trial.model.label())

        if add_scatter:
            self.ax[0].scatter(*trial.xy("ZnO", pct=True))

        self.ax[1].plot(*trial.model.rr_v_y())
        self.ax[1].set_xlim(12, 0)
        self.ax[1].set_ylim(0.01, 10)
        self.ax[1].set_yscale("log")
        self.ax[1].set_xlabel("ZnO [mol-%]")


class DetailPlot(Plot):

    def __init__(self, setup, trial_id, excludes=[],
                 y0=None, title="", line_styles=["solid", "dashed", "dotted"],
                 colors="k", markers=["o", "s"],
                 t=None, y_max=12, axes=[], show_ylabel=False):
        super(DetailPlot, self).__init__(setup)
        trial = SlagReductionTrial(VxPySelector(trial_id)(self.phased_data))
        self.t = t
        self.trial = self.assembler.do(trial, exclude=excludes)
        self.exclude = excludes
        self.y0 = y0

        self.line_styles = line_styles
        self.colors = colors
        self.markers = markers
        self.cnt = 0

        if len(axes) == 0:
            self.fig, self.ax = plt.subplots(nrows=2)
            self.fig.subplots_adjust(wspace=0, hspace=0)
        elif len(axes) == 2:
            self.ax = axes
        else:
            raise ValueError("Length of axes must be two.")

        #self.ax_model.set_xlabel("Time [min]")

        self.ax_model.set_title(title)
        self.ax_model.set_ylim(0, y_max)
        self.ax_dev.plot([min(t), max(t)], [0, 0], color="lightgrey")
        self.ax_dev.set_ylim(-0.99, 0.99)
        self.ax_dev.set_xlabel("Time [min]")

        if show_ylabel is True:
            self.ax_model.set_ylabel("ZnO [mol-%]")
            self.ax_dev.set_ylabel("$x\displaystyle_{ZnO,Model}-x_{ZnO,Experimental}$ [mol-%]")

        else:
            self.ax_model.set_yticks([])
            self.ax_dev.set_yticks([])

    @property
    def ax_model(self):
        return self.ax[0]

    @property
    def ax_dev(self):
        return self.ax[1]

    def scatter(self, which):
        compound = self.trial.model_compounds.name_x(which)
        self.ax_model.scatter(*self.trial.xy(compound, pct=True),
                              marker=self.markers[self.cnt],
                              c=self.colors)

    def add(self, model, which=0, color=None, line_styles=None):
        t = self.t
        self.assembler.model = model
        self.trial.model = model
        trial = self.assembler.do(copy.copy(self.trial))

        if color is not None:
            self.colors = color

        if line_styles is not None:
            self.line_styles = line_styles

        if self.y0 is not None:
            self.trial.model.shift_to_y0(self.y0)

        if type(self.colors) is str:
            color = self.colors
        elif type(self.colors) is list:
            color = self.colors[self.cnt]
        else:
            raise TypeError("Colors must be a list or str")

        self.ax_model.plot(*trial.model_y(pct=True, t=t, which=which),
                               label=trial.model.label(which),
                               linestyle=self.line_styles[self.cnt], color=color)

        self.ax_dev.plot(*trial.model.deviation(which=which),
                         linestyle=self.line_styles[self.cnt],
                         color=color, label=f"{trial.model.fit_report}")

        self.cnt += 1

    def legend(self):
        self.ax_model.legend(fontsize="small")
        self.ax_dev.legend(fontsize="small")
