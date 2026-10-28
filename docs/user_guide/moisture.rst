水汽：相对湿度、露点和可降水量
==============================

相对湿度和露点用的是 UPP 自己的饱和公式，不是 Bolton，也不是 Magnus 的常见写法。常数在 ``params.F``：``PQ0 = 379.90516``，``A2 = 17.2693882``，``A3 = 273.16``，``A4 = 35.86``。

相对湿度
--------

``CALRH.f``：

.. math::

   q_s = \frac{P_{q0}}{p}
     \exp\left(A_2\frac{T - A_3}{T - A_4}\right), \quad
   \mathrm{RH} = q / q_s

结果限制在 ``[10^{-6}, 1]``。碰到上下界时，Fortran 会把比湿改成 ``\mathrm{RH}\,q_s``。气压层插值里这一步是打开的；2 m 相对湿度先算分数再乘 100 写成百分数，和 GRIB 一致。

这支是 WRF 用的 ``CALRH``，不是 GFS 的 ``CALRH_GFS``，也不是 RAP 的 ``CALRH_GSD``。

露点
----

``CALDWP.f`` 先算水汽压（帕），再换成厘巴：

.. math::

   e = \frac{p q}{\epsilon + (1-\epsilon) q}, \quad
   \epsilon = 18.015 / 28.964

``DEWPOINT.f`` 用下面的式子造一张 2000 点的表，范围 0.0001 到 10 厘巴，再用线性插值：

.. math::

   e = 0.611\, X^{A} \exp\big((A+B)(1-X)\big), \quad
   X = 273.16 / T_d

``A = 5.0065``，``B = 19.83923``。表用牛顿法逐点反演，和 Fortran 的循环相同。最后把露点限制成不超过环境温度。

可降水量
--------

``CALPW.f`` 在 ``IDECID <= 1`` 时只积分水汽，不加水凝物：

.. math::

   \mathrm{PWAT} = \sum_k q_k\,\Delta p_k / g

``\Delta p`` 是界面气压差。单位 kg m\ :sup:`−2`，和毫米水当量相同。云水、云冰有单独的 ``IDECID``，这里没有当作 PWAT。

.. code-block:: python

   from upp_py.thermo import dewpoint, relative_humidity
   from upp_py.post import precipitable_water

   rh = relative_humidity(state.p_shelter, state.t2, state.q2, percent=True)
   td = dewpoint(state.p_shelter, state.q2, state.t2)
   pwat = precipitable_water(state)

和 UPP 的差别
-------------

CAPE 里的饱和水汽压用的是另一套，``fpvsnew`` 的 Clausius–Clapeyron（水相 273.16 K 以上，冰相 253.16 K 以下，中间线性混合）。那是 ``CALCAPE.f`` 算浮力时的选择，相对湿度和露点不要误用它。
