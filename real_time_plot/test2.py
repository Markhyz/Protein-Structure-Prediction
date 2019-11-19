# import matplotlib.pyplot as plt

# plt.plot([1, 2], [2, 3])

# plt.show()

import threading

def inp():
  while 1:
    x = input()
    print(x)

t = threading.Thread(target=inp)

t.start()

while 1:
  for i in range(100000000):
    continue
  print("lol")