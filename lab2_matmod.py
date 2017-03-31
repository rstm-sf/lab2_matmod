import numpy as np
import matplotlib.pyplot as plt

def F(N, kd, k, dt, l, l0, m, D, nu, ki):
	x = np.empty(N); y = np.empty(N);

	x[0] = l0;  y[0] = 0.0;
	x[N-1] = l; y[N-1] = 0.0;

	s = 0.0; dt_m = dt/m; two_D = 2.0/D;
	for i in range(N-1):
		x_i = x[i]; y_i = y[i];
		x[i+1] = y_i * dt + x_i;
		s   += (2 * l - x[i+1] - x_i) * dt * 0.5;
		y[i+1] = dt_m * (two_D * ( (l - x_i)*k + kd*y_i + ki*s ) - nu*y_i) + y_i;

	return x, y


l = 1.0; l0 = 0.2; m = 1.0; D = 6.0; nu = 0.5;
kd1 = 0; k1 = 0.6; ki1 = 0;
N = 400; dt = 0.05;
t = [i*dt for i in range(N)]

fig  = plt.figure(figsize=(10, 8), dpi=100)
subplot1 = fig.add_subplot(221)   #left
subplot2 = fig.add_subplot(222)   #right

xmin = float("inf"); ymin = float("inf");
xmax = -float("inf"); ymax = -float("inf");

for pp in range(4):
	ki = ki1 - 0.01*pp

	for i in range(1):
		kd = kd1 - 0.5*i

		for j in range(1):
			k = k1 - 0.1*j
			x, y = F(N, kd, k, dt, l, l0, m, D, nu, ki)

	subplot1.plot(t, x)
	subplot2.plot(t, y)

	if xmin > x.min():
		xmin = x.min()
	if ymin > y.min():
		ymin = y.min()
	if xmax < x.max():
		xmax = x.max()
	if ymax < y.max():
		ymax = y.max()

subplot1.set_ylim(xmin, xmax)
subplot2.set_ylim(ymin, ymax)
subplot1.grid()
subplot2.grid()
plt.show()