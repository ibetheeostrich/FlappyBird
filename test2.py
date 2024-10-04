import numpy as np
import matplotlib.pyplot as plt
import matplotlib

from aero_solver.bem import bem_span as bem_span

t_step = 0.00125
no_steps = 400
chords = np.zeros(no_steps)#*np.linspace(0,0.25,no_steps)
freq = 3
amp = np.deg2rad(42.5)


tag = 1

_, cl, cd ,td = bem_span(tag, chords, t_step, no_steps, freq, amp)

# Creating Figures
fig, ax = plt.subplots()
fig.dpi = 300
ax.plot(td, cl)
ax.plot(td, cd)
ax.set_xlabel('Time  (s)')
ax.set_ylabel('Lift Force (n)')
plt.savefig('long_test' + '.png')
plt.close(fig)
