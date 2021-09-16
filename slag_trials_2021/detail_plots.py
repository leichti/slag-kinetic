from nemlib.model import ZnOFeModel, FirstOrderModel, ZnOCCrucibleModel, NickelModel
from setup import MetalBathSeries, SummaryPlot, DetailPlot
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
matplotlib.rcParams['text.usetex'] = True

fig, axes = plt.subplots(figsize=(12,7), ncols=4, nrows=2)
fig.subplots_adjust(hspace=0, wspace=0, left=0.08, right=0.98)

plot = DetailPlot(MetalBathSeries(), "V27", [0, 1],
                  title="Graphite Crucible",
                  t=np.arange(0, 210, 1), axes=axes[0:, 0],
                  show_ylabel=True)
plot.scatter(0)
plot.add(model=ZnOCCrucibleModel)
plot.add(model=FirstOrderModel)
plot.legend()

plot = DetailPlot(MetalBathSeries(), "V34", [0, 1],
                  title="+ Copper",
                  t=np.arange(0, 200, 1), axes=axes[0:, 1])
plot.scatter(0)
plot.add(model=ZnOCCrucibleModel)
plot.add(model=FirstOrderModel)
plot.legend()

plot = DetailPlot(MetalBathSeries(), "V28", [0, 1],
                  title="+ Nickel",
                  t=np.arange(0, 75, 1), axes=axes[0:, 2])
plot.scatter(0)
plot.add(model=NickelModel, )
plot.add(model=ZnOCCrucibleModel)
plot.add(model=FirstOrderModel)
plot.legend()

plot = DetailPlot(MetalBathSeries(), "V23", [0],
                  title="+ Iron",
                  t=np.arange(0, 110, 1), axes=axes[0:, 3])
plot.scatter(0)
plot.add(model=ZnOFeModel, which=0)
plot.scatter(1)
plot.add(model=ZnOFeModel, which=1)
plot.legend()

fig.savefig("plots/detail_plots.svg")
plt.show()