import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import threading

x = []
y = []

def inp():
  while 1:
    valx, valy = [int(val) for val in input().split()]
    x.append(valx)
    y.append(valy)
    with open("lel", "w") as file:
      file.write("{} {}\n".format(x, y))

t = threading.Thread(target=inp)
t.start()

figure = plt.gcf()
line, = plt.plot([], [])

def update(frame):
  line.set_data(x, y)
  if len(x) > 0:
    figure.gca().relim()
    figure.gca().autoscale_view()
  return line,

animation = FuncAnimation(figure, update, interval=100)

plt.show()