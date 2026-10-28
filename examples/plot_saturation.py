"""
饱和比湿与露点
================

``CALRH.f`` 的饱和比湿，以及 ``DEWPOINT.f`` 在 1000 hPa 上对应的露点。
这是气压层相对湿度和 2 m 露点共用的公式。
"""

import matplotlib.pyplot as plt
import numpy as np

from upp_py.thermo import dewpoint, saturation_specific_humidity

t = np.linspace(250.0, 310.0, 200)
p = np.full_like(t, 1.0e5)
qs = saturation_specific_humidity(t, p)
# 取饱和比湿的一半，看露点落在哪里。
q = 0.5 * qs
td = dewpoint(p, q, t)

fig, ax = plt.subplots(figsize=(6.2, 3.8))
ax.plot(t, qs * 1000.0, label=r"$q_s$ (g/kg)")
ax.set_xlabel("Temperature (K)")
ax.set_ylabel("Specific humidity (g/kg)")
ax.set_title("UPP saturation specific humidity at 1000 hPa")
ax.legend()
fig.tight_layout()
plt.show()

fig, ax = plt.subplots(figsize=(6.2, 3.8))
ax.plot(t, td, label="dewpoint at q = 0.5 qs")
ax.plot(t, t, color="0.6", linewidth=0.8, label="T = Td")
ax.set_xlabel("Temperature (K)")
ax.set_ylabel("Dewpoint (K)")
ax.legend()
fig.tight_layout()
plt.show()
