能见度
======

对应 ``CALVIS.f``。``MDLFLD.f`` 在最低质量层调用它，可用场 180，短名 VIS，单位米。

晴空时消光只剩 ``1e-10`` km\ :sup:`−1`。``-ln(0.02)/beta`` 会非常大，写出前封在 24.135 km，也就是 24135 m。有水成物时四项相加：

- 云水：``144.7 * C**0.88``
- 雨水：``2.24 * C**0.75``
- 云冰：``327.8 * C``
- 雪：``10.36 * C**0.7776``

``C`` 是 g m\ :sup:`−3`。系数来自 Kunkel、Marshall-Palmer、Rutledge 与 Hobbs、Stallabrass，和 Fortran 注释一致。

水汽用 ``INITPOST.F`` 转成比湿之后的 ``Q``。``CALVIS`` 的注释写混合比，调用处传的是这个 ``Q``。

``RAINCV > 0`` 时，先按 ``MDLFLD.f`` 把对流降水折进最低层的雨水混合比，再算能见度。``RAINCON = 0.8333*1.1787e4``，``SNOCON = 0.94*1.4594e5``。没有 ``SR`` 时这一步全部算成雨。``imp_physics == 99`` 是 NMMB 的 Zhao 方案，WRF-ARW 不走那一支。

``CALVIS_GSD``（可用场 410）只在 ``imp_physics == 28`` 时调用，还要 ``ZENSUN.f`` 的太阳天顶角余弦。这里不算。短名 ``VIS_GSD`` 会写进数据集属性 ``skipped``。
