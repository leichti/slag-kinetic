from nemlib.model import FirstOrderModel
from setup import plot_trial
import matplotlib.pyplot as plt

plot_trial("V34", [0, 1], model=FirstOrderModel, add_scatter=True)
#plot_trial("V28", [0, 1])
#plot_trial("V27", [0, 1])
#plot_trial("V23", [0])

plt.legend()
plt.show()

