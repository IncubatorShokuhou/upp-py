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

桶和多时次
----------

``SURFCE.f`` 的注释假定累积量不会把桶倒空，写出的 APCP 就是 ``RAINC+RAINNC``。WRF 在 ``bucket_mm > 0`` 时会把 ``RAINC``/``RAINNC`` 减去 ``bucket_mm``，并给 ``I_RAINC``/``I_RAINNC`` 加 1。文件里有计数器和正的 ``BUCKET_MM`` 时，连续总量是

.. math::

   \mathrm{APCP} = \mathrm{RAINC}+\mathrm{RAINNC}+B\,(I_{c}+I_{nc})

``B`` 是 ``bucket_mm``，``I_c`` 和 ``I_nc`` 是 ``I_RAINC`` 和 ``I_RAINNC``。

没有计数器时仍用 ``RAINC+RAINNC``，和 ``unipost`` 写出的数相同。

另一套桶是 PREC_ACC_C 和 PREC_ACC_NC，单位毫米，时间长度是全局属性 PREC_ACC_DT，单位秒。INITPOST.F 把它们加进 PCP_BUCKET。SURFCE.f 可用场 434 在预报时效为 0 时写成 0，其它时效写成这个和。对流桶是可用场 435，格点尺度桶是 436。雪桶 SNOW_ACC_NC（可用场 437，GRIB 名 WEASD）在时效 0 不清零。

短名 APCP 是从起报开始的连续总量。APCP_BUCKET 是上面的 PREC_ACC 桶。命令行 ``accum-hours`` 再写变量 apcp_accum：窗口终点是每个选中的时效，起点早于起报时收到起报时刻，和 SURFCE.f 把时间范围起点收到 0 一样。历史时次对不齐窗口起点时，用不晚于起点的最近一个时次，实际小时数在坐标 accum_hours_actual。

没有 ``I_RAINC`` 而 ``RAINC`` 在两个时次之间下降时，下降的格点改用后一时次的 ``PREC_ACC`` 桶，当作这一段的增量。

.. code-block:: bash

   upp-py wrfout_d01 -o post.nc --times 0,1,2,3 --accum-hours 1

和 UPP 的差别
-------------

``unipost`` 一次只处理 itag 里的一个时次，并且不把 ``I_RAINC*bucket_mm`` 加回去。多时次窗口和计数器修正是在这份 wrfout 上把桶倒空补上，不是 ``SURFCE.f`` 里现成的循环。参考 GRIB 的起报时刻 APCP 仍是 0。
