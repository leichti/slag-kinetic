from nemlib.model import ZnOFeModel, FirstOrderModel, ZnOCCrucibleModel, NickelModel
from setup import MetalBathSeries, SummaryPlot
import matplotlib.pyplot as plt


plot = SummaryPlot(MetalBathSeries())
plot.add("V27", [0,1], model=ZnOCCrucibleModel, add_scatter=True)
plot.add("V34", [0,1], model=ZnOCCrucibleModel, add_scatter=True)
plot.add("V28", [0,1], model=NickelModel, add_scatter=True)
plot.add("V23", [0], model=ZnOFeModel, add_scatter=True)

plot.ax[0].legend()

plt.show()

