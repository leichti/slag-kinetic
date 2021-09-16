from nemlib.model import ZnOFeModel, FirstOrderModel, ZnOCCrucibleModel
from setup import MetalBathSeries, SummaryPlot
import matplotlib.pyplot as plt

plot = SummaryPlot(MetalBathSeries())
plot.add("V23", [0], model=ZnOFeModel, add_scatter=True)


plt.legend()
plt.show()

