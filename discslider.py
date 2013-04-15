import matplotlib.pyplot as plt
import numpy as np
from matplotlib.widgets import Slider

class ChangingPlot(object):
    def __init__(self):
        self.inc = 0.5

        self.fig, self.ax = plt.subplots()
        self.sliderax = self.fig.add_axes([0.2, 0.02, 0.6, 0.03],
                                          axisbg='yellow')

        self.slider = DiscreteSlider(self.sliderax, 'Value', 0, 10, 
                                     increment=self.inc, valinit=self.inc)
        self.slider.on_changed(self.update)

        x = np.arange(0, 10.5, self.inc)
        self.ax.plot(x, x, 'ro')
        self.dot, = self.ax.plot(self.inc, self.inc, 'bo', markersize=18)

    def update(self, value):
        self.dot.set_data([[value],[value]])

    def show(self):
        plt.show()

class DiscreteSlider(Slider):
    """A matplotlib slider widget with discrete steps."""
    def __init__(self, *args, **kwargs):
        """Identical to Slider.__init__, except for the "increment" kwarg.
        "increment" specifies the step size that the slider will be discritized
        to."""
        self.inc = kwargs.pop('increment', 0.5)
        Slider.__init__(self, *args, **kwargs)

    def set_val(self, val):
        discrete_val = float(val / self.inc) * self.inc
        # We can't just call Slider.set_val(self, discrete_val), because this 
        # will prevent the slider from updating properly (it will get stuck at
        # the first step and not "slide"). Instead, we'll keep track of the
        # the continuous value as self.val and pass in the discrete value to
        # everything else.
        xy = self.poly.xy
        xy[2] = discrete_val, 1
        xy[3] = discrete_val, 0
        self.poly.xy = xy
        self.valtext.set_text(self.valfmt % discrete_val)
        if self.drawon: 
            self.ax.figure.canvas.draw()
        self.val = val
        if not self.eventson: 
            return
        for cid, func in self.observers.iteritems():
            func(discrete_val)


#p = ChangingPlot()
#p.show()
