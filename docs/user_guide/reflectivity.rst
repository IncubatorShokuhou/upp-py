反射率
======

对应 ``MDLFLD.f`` 里 Thompson（``mp_physics = 8``）那一支，以及它前面调用的 ``PARAMR.f``。其他方案走同一文件里的 ``(q\rho)^{1.75}`` 估算。

Thompson / Stoelinga
--------------------

雨、雪、霰用指数谱。雨的反射率因子（先以 m\ :sup:`6` m\ :sup:`−3` 计，再乘 ``10^{18}`` 换成 mm\ :sup:`6` m\ :sup:`−3`）：

.. math::

   Z_r = 720\, N_{0r}\, \lambda_r^{-7}

截距 ``N_{0r}`` 在毛毛雨和常规雨之间用双曲正切过渡（``PARAMR.f`` 的 ``qr0 = 2\times 10^{-4}``）。雪的截距随温度变化，取 Houze 等 1979 的形式，上限 ``2\times 10^8`` m\ :sup:`−4`。霰用 Rasmussen 的变截距。冰相再乘 ``\alpha = 0.224`` 和密度比的平方。亮带放大系数 ``bb`` 在这版源码里是 0。

对流贡献是 ``Z = 300 R^{1.4}``，``R`` 为毫米每小时，并且不小于 0.1 mm\ :sup:`6` m\ :sup:`−3`。所以没有水凝物时诊断值大约是 −10 dBZ。

写出的合成反射率
----------------

NCAR（以及 NMM、FV3）加上 Thompson 时，写进 GRIB 的合成反射率不是上面这个诊断，而是模式变量 ``REFL_10CM`` 的柱最大，再限制在 ``[-20, 80]`` dBZ。Carr Fire 的 wrfout 里 ``REFL_10CM`` 全是 0，参考 GRIB 的 ``REFC`` 因此也是 0。unipost 日志里的 ``dbze_max`` 是诊断值，大约 1 dBZ，没有写到这条记录上。

1 km 和 4 km 高度上，非 RAP 的分支用的是诊断值：从地面往上，第一个顶界面超过该高度的质量层。RAP/HRRR 才会直接读 ``REFL_10CM_1KM``。这个个例的参考 GRIB 把 1 km 反射率写成了 0，和诊断值的 −10 dBZ 对不上，验证页里单独说明了。

其他微物理
----------

``mp_physics`` 不是 8 时，用

.. math::

   Z = (q\rho)^{1.75} \times
   \begin{cases}
   3.630803\times 10^{-9} & \text{雨} \\
   2.18500\times 10^{-10} & \text{雪} \\
   1.033267\times 10^{-9} & \text{霰}
   \end{cases}
   \times 10^{18}

密度 ``\rho = p / (R_d T (1 + 0.608 q))``。

.. code-block:: python

   from upp_py.reflectivity import reflectivity

   out = reflectivity(state)
   out.composite   # 按 NCAR 分支写出的 REFC
   out.diagnosed   # Stoelinga 柱最大
   out.agl_1km

和 UPP 的差别
-------------

诊断值和查表无关，公式是直接算的。模式自己的 ``REFL_10CM`` 若是空的，合成反射率会像参考场一样是 0；需要诊断回波时看 ``diagnosed``，不要看写出的 ``REFC``。
