降水
====

对应 ``INITPOST.F`` 读累积量、``SURFCE.f`` 写出的那一段。

累积降水
--------

WRF 的 ``RAINC`` 和 ``RAINNC`` 单位是毫米。UPP 先乘 0.001 变成米，写 GRIB 时再乘 1000，回到 kg m\ :sup:`−2`。所以

.. math::

   \mathrm{APCP} = \mathrm{RAINC} + \mathrm{RAINNC}

桶降水 ``PREC_ACC_C`` 和 ``PREC_ACC_NC`` 在预报时效为 0 时被置 0（``SURFCE.f`` 里 ``IFHR == 0`` 的分支）。这个个例的有效时刻就是起报时刻，累积量本身也是 0。

降水率
------

瞬时总降水率用当前步的 ``RAINCV + RAINNCV``（毫米）和模式步长 ``DT``：

.. math::

   \mathrm{PRATE} = (\mathrm{RAINCV}+\mathrm{RAINNCV})\times 10^{-3} / \Delta t \times 1000

也就是毫米每秒，和 kg m\ :sup:`−2` s\ :sup:`−1` 相同。对流率单独用 ``RAINCV``。Thompson 反射率里的对流降水率还要再乘 3600 换成毫米每小时，那是 ``MDLFLD.f`` 的 ``CUPRATE = 3.6\times 10^6/\Delta t \times \mathrm{CPRATE}``，``CPRATE`` 已经是米每步。

.. code-block:: python

   apcp = state.rainc + state.rainnc          # kg/m2
   prate = (state.raincv + state.rainncv) * 0.001 / state.dt * 1000.0

和 UPP 的差别
-------------

多时次之间的桶清空、3 小时或 6 小时累积窗口没有做。只有一个时次时，APCP 就是从积分开始到该时次的总量。参考 GRIB 这个时次的 APCP 是 0，两边一致。
