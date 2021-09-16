from nemlib.model import ZnOFeModel, FirstOrderModel, ZnOCCrucibleModel
from setup import MetalBathSeries, SummaryPlot, DetailPlot
import numpy as np

plot = DetailPlot(MetalBathSeries(), "V27", [0, 1], title="Carbon Crucible")
plot.add(model=ZnOCCrucibleModel, t=np.arange(0,250,1))
plot.add(model=FirstOrderModel, t=np.arange(0,250,1))

plot.finalize()


