气压层与高度插值
================

对应 ``MDL2P.f``。标准气压层是 ``CTLBLK.f`` 的 ``SPLDEF``，从 2 hPa 到 1000 hPa，一共 46 层。配置里的 ``pressure_hpa`` 可以是它的子集。膜法海平面气压仍然会在内部把这 46 层都插出来，因为积分要从上往下走。

对数气压插值
------------

先找目标气压正下方的质量层 ``L``（气压更大的那一层）和它上面的 ``L-1``：

.. math::

   f = f_L + (f_L - f_{L-1})
     \frac{\ln p_t - \ln p_L}{\ln p_L - \ln p_{L-1}}

NCAR 分支把这个比例限制在 ``[-10, 10]``。目标层正好落在某个质量层上时，比例是 0 或 −1，结果就是该层的值。

三种垂直位置：

1. 目标气压高于最上面两层之间的界面：温度、比湿和风直接抄模式顶。
2. 最低质量层以下、地面以上：仍用上面的公式。比例为正，沿最底下两层的斜率往下延。地面气压只比目标层高几帕时走的就是这一支，温度仍接近最低质量层。
3. 目标气压大于等于地面气压：走地下外推。每个格点都写出数，等压面图上不留空。

地面以上的比湿会按 ``CALRH`` 的饱和比湿把相对湿度夹到 ``[10^{-6}, 1]``，再写回比湿。

地下外推
--------

参考状态取离地面第二层和第三层的平均，气压用这两层之间的界面 ``PINT``（WRF 顺序里是 ``p_int[2]``）。相对湿度先夹到 ``[10^{-6}, 1]``，外推时保持这个相对湿度。虚拟温度按 ``params.F`` 的标准大气递减率走：

.. math::

   \gamma = 6.5\times 10^{-3}\,\mathrm{K\,m^{-1}}, \quad
   R_d\gamma/g = \mathtt{RGAMOG}

   T_v(p_t) = T_v(p_L)\,(p_t / p_L)^{R_d \gamma / g}

比湿用外推后的温度重算，下限 ``10^{-12}``。风和垂直速度抄最低质量层。水凝物在地下置零。

地下的位势从上一张标准气压面往下做一步静力积分，所以气压层必须从低压到高压处理：

.. math::

   \Phi(p) = \Phi(p_\mathrm{prev})
     - \frac{R_d}{2}
       \big(T_v(p_\mathrm{prev}) + T_v(p)\big)
       \ln\frac{p}{p_\mathrm{prev}}

上一张气压面若还在地面以上，它的高度来自界面高度的 ``\ln p`` 插值，再乘 ``g`` 得到位势。这样跨过地面的那一步只跨一个标准层距（通常 25 hPa），后面的地下层再一层层接下去。膜法海平面气压内部会把 ``SPLDEF`` 的 46 层都算出来，输出子集也沿用这套高度。

地上的高度不走这条积分。非 GFS 的分支在界面高度 ``(PH+PHB)/g`` 上做同样的 ``\ln p`` 插值，高度上的比例不裁剪。界面气压除了地面和模式顶，是相邻两个质量层气压的平均。判断地下用的地面气压是 ``INITPOST.F`` 里加湿、平滑之后的 ``PINT(LM+1)``。

绝对涡度在插到气压层之后的风上计算，对应 ``MDL2P.f`` 里对 ``USL, VSL`` 调用 ``CALVOR``。A 网格的公式是

.. math::

   \zeta_a = \frac{\partial v}{\partial x} - \frac{\partial u}{\partial y}
     + f + \frac{u\tan\phi}{a}

``\Delta x = DX / m``，``a = 6.371\times 10^6\,\mathrm{m}``。最外一圈 Fortran 留成初值 0，这里用内侧一点填上，避免图上出现一圈零。

.. code-block:: python

   from upp_py.interp import interpolate_pressure_levels
   from upp_py.wrf import open_wrfout

   state = open_wrfout("wrfout_d01")
   levels = interpolate_pressure_levels(state, [85000.0, 50000.0, 25000.0])
   levels[1].temperature  # 500 hPa，返回顺序是自上而下

和 UPP 的差别
-------------

GFS 在地下用 Shuell 递减率，WRF-ARW 走上面的 NCAR 分支。850 hPa 在高山上会落到地下，外推之后每个格点都有高度、温度和湿度。
