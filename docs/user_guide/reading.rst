读取 wrfout 与去交错
=======================

对应 ``INITPOST.F`` 和 ``getVariable.f``。WRF-ARW 在 UPP 里的模式名是 ``NCAR``，网格类型是 A 网格：质量点上的风是 U、V 去交错之后的结果。

垂直方向
--------

WRF 的 ``bottom_top`` 从地面往上数。UPP 读进来之后把 k 翻过来，让 ``L = 1`` 在模式顶、``L = LM`` 在最低层，因为后面的气压插值假定下标增大时气压升高。upp-py 在内存里保持 WRF 的顺序（``k = 0`` 在地面），公式按翻转后的下标改写，避免两套数组来回倒。

``getVariable.f`` 里那句 ``ll = lm1 - l + 1`` 就是这次翻转。

去交错
------

质量点上的风和垂直速度（``INITPOST.F``）：

.. math::

   u(i) &= \tfrac12 \big(U(i) + U(i+1)\big) \\
   v(j) &= \tfrac12 \big(V(j) + V(j+1)\big) \\
   w(k) &= \tfrac12 \big(W(k) + W(k+1)\big)

U 在 x 方向多一个点，V 在 y 方向多一个点，W 和位势在垂直方向多一个点。

热力学状态
----------

位温扰动存在 ``T`` 里。UPP 先加 300 K，再用气压把它换成温度：

.. math::

   T = (\theta' + 300)\,(p / 10^5)^\kappa, \quad \kappa = 0.28589641

``\kappa`` 是 ``params.F`` 里的 ``CAPA``，和 ``R_d/c_p`` 差一点点，不要换成后者。

v4.1 的气压用 ``P_HYD``，不再把 ``P`` 和 ``PB`` 加起来。文件里没有 ``P_HYD`` 时才退回 ``P+PB``。比湿由混合比换算，``q = r / (1+r)``。水凝物缺测时按 0，和 ``getVariable`` 找不到变量就填 0 一样。

位势高度是 ``(PH + PHB) / g``，``g = 9.81``。质量层高度是上下两个界面的平均，也就是 ``ZMID``。

地面气压不是直接拿 ``PSFC``。非 RAP 的分支用干空气柱 ``MU + MUB + P_TOP``，再加上整层水汽换算成的气压，并用 8 点平滑迭代 405 次（Pyle/Chuang 对 ARW 的修正）。平滑系数来自 ``OVERRC = 1.10``。

投影
----

从全局属性读 ``MAP_PROJ``、``TRUELAT1`` 和 ``TRUELAT2``，以及 ``STAND_LON``、``CEN_LAT``，还有 ``DX`` 与 ``DY``。1 是 Lambert，2 是极射赤面，3 是 Mercator。``MAPFAC_M`` 用来把名义格距换成当地格距：``\Delta x = DX / m``。经纬度用质量点上的 ``XLAT``、``XLONG``。

有效时间来自 ``Times``，起报时间来自 ``START_DATE``。预报时效是两者之差，这个个例的起始时刻和有效时刻相同，所以是 fhr 00。

.. code-block:: python

   from upp_py.wrf import open_wrfout

   state = open_wrfout("wrfout_d01_2018-07-23_21_00_00")
   state.t.shape   # (nz, ny, nx)，k=0 在地面
   state.p_sfc     # 含平滑水汽订正的地面气压

和 UPP 的差别
-------------

UPP 用 WRF I/O API 按记录读，这里用 xarray 一次读完。缺测变量填 0 的行为和 ``getVariable`` 相同。土壤层、气溶胶和 GOES 亮温不在这个版本的范围内。
