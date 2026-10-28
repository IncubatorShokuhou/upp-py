对流指数：CAPE 与 CIN、LCL 和 LFC
==========================================

对应 ``CALCAPE.f``，气块的选法在 ``CLDRAD.f``、``BNDLYR.f`` 和 ``MISCLN.f``。

积分
----

从抬升凝结高度积到平衡高度（气块比环境暖的最高一层）。正的面积进 CAPE，负的进 CIN，CIN 保持为负数：

.. math::

   \sum g \ln\frac{\theta_{v,p}}{\theta_{v,a}}\,\Delta z

虚位温用 ``fpvsnew`` 的饱和水汽压算气块比湿，环境用实际比湿。``\kappa`` 仍是 ``CAPA``。没有正面积时 CAPE 和 CIN 都是 0，因为平衡层会留在地面，求和区间是空的。

气块相当位温按 ``CALCAPE.f`` 写成

.. math::

   \theta_e = \theta \exp\left(
     \frac{L}{c_p}\, q\, (10^5/p_{LCL})^\kappa / \theta
   \right)

``L/c_p`` 是 ``ELOCP = 2.72\times 10^6 / 1004.6``。``p_{LCL}`` 是干绝热线上饱和比湿等于气块比湿的气压。UPP 用 ``TABLE.f`` 的 76×134 表来查，这里对同一公式做二分。湿绝热上的温度是把这个 ``\theta_e`` 反解出来，对应 ``TTBLEX``，同样不经过样条表。直接求解是有意的：优先方程本身，不把查表量化对齐到 ``unipost.exe``。

四种气块
--------

.. list-table::
   :header-rows: 1

   * - 名字
     - GRIB 层
     - 选法
   * - ``surface``
     - surface
     - ``CLDRAD.f``：最低 10 hPa 里 ``\theta_e`` 最大的模式层。注释写 70 hPa，代码是 ``DPBND = 10e2``
   * - ``mixed_layer``
     - 90–0 hPa
     - 地面往上三个 30 hPa 层的平均，``BNDLYR.f``
   * - ``best``
     - 180–0 hPa
     - 六个 30 hPa 层里 Bolton ``\theta_e``（``CALTHTE.f``）最大的一层
   * - ``most_unstable``
     - 标成 255–0 hPa
     - ``MISCLN.f``：地面往上 300 hPa 里 ``\theta_e`` 最大的模式层

边界层的层平均是气压厚度加权，层的代表气压是上下界面的平均，不是质量加权气压。

LCL 是饱和气压，高度在界面高度上按 ``\ln p`` 插到这个气压，再减地形。LFC 是 LCL 之上第一层浮力变正的质量层。

.. code-block:: python

   from upp_py.cape import cape_mixed_layer, cape_surface

   sfc = cape_surface(state)
   sfc.cape, sfc.cin, sfc.lcl_height_agl

和 UPP 的差别
-------------

upp-py 对湿绝热和 Bolton :math:`\theta_e` 直接求解，有意不用 ``CALCAPE`` / ``TTBLEX`` 的查表。和 ``unipost.exe`` 比，最佳 CAPE 剩下的差是预期结果。六个 30 hPa 层里 Bolton :math:`\theta_e` 非常接近时，两边可能选到相邻的一层，CAPE 会差一截。Carr Fire 上这种情况大约占 0.2% 的格点，所以“最佳 CAPE”的 RMSE 比地面 CAPE 大，99 分位绝对差仍然在 20 J/kg 附近。地面、混合层和最不稳定 CAPE 的 RMSE 大约是 10 J/kg。
