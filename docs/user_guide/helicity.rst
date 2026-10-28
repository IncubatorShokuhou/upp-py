螺旋度与风暴移动
================

对应 ``CALHEL.f``。风暴移动用 Bunkers 等 1998 的方法，不再用 Davies-Jones 1993 那种把平均风偏转 20° 或 30° 的做法。

0–6 km、0–0.5 km 和 5.5–6 km 的平均风是质量点风速的算术平均，不按质量加权。高度是离地高度，地形用 ``HGT``（``FIS/g``）。5.5–6 km 没有层时，改拿 6–7 km 里从地面往上遇到的第一层。

风暴移动在 0–6 km 平均风右侧 7.5 m/s，方向垂直于（5.5–6 km 风减去 0–0.5 km 风）：

.. math::

   \mathbf{c} = \mathbf{v}_{0-6}
     + 7.5\,\frac{(-v_s,\, u_s)}{|\mathbf{s}|}

``\mathbf{s}`` 是上面的切变。A 网格上四个相邻点的平均退回成该点自己，WRF-ARW 就是这种网格。

螺旋度在 0–1 km 和 0–3 km 两段里，用上下层风速差估计垂直切变，再和相对风暴的风做叉乘，乘上层厚度。写法与 ``CALHEL.f`` 里 ``HELI`` 的累加相同。

.. code-block:: python

   from upp_py.helicity import storm_relative_helicity

   h = storm_relative_helicity(state)
   h.hel_0_1km, h.hel_0_3km
   h.u_storm, h.v_storm

和 UPP 的差别
-------------

Carr Fire 上两个深度的螺旋度和参考 GRIB 的差在 0.05 m\ :sup:`2` s\ :sup:`−2` 以内，就是打包误差。边缘点如果三层平均里有一层数不到风，风暴移动按 Fortran 置 0。
