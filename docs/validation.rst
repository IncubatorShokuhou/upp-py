与 unipost.exe 的对比
======================

个例是 Carr Fire，``wrfout_d01_2018-07-23_21_00_00``，WRF-ARW 4.7.1，199×199×79，3 km Lambert，有效时间 2018-07-23 21:00 UTC，和起报时间相同。参考场是 UPP v4.1 ``unipost.exe``（gfortran）写出的 ``WRFPRS_d01.00``，549 条 GRIB2。

差是 upp-py 减去参考场。重新生成表格：

.. code-block:: bash

   python scripts/validate_against_upp.py \
       --wrfout /path/to/wrfout_d01_2018-07-23_21_00_00 \
       --grib /path/to/WRFPRS_d01.00

同一份数字也写在 ``docs/validation.md``。

.. include:: validation_generated.rst

怎么读这张表
------------

PRMSL 是 ``NGMSLP.f`` 的 Shuell 订正（可用场 105）。膜法是 MSLET，``SLP_new.f`` 的 ``MEMSLP``（可用场 23）。日志里的 ``CALLING MEMSLP`` 对应的是 MSLET，不是 PRMSL。两条都和各自的 GRIB 记录比过，均方根差大约 6 Pa。

2 m 的温度、露点和相对湿度，加上 10 m 风，以及 250 hPa、500 hPa 的高度和温度，差在 GRIB 打包的量级。850 hPa 在高山上有一部分在地下：温度用第二、三层的平均虚温按 ``RGAMOG`` 外推，高度从上一张标准气压面做一步静力积分，每个格点都有值。高度最大差约 0.01 m，温度约 0.05 K，相对湿度约 0.5%。

地面、混合层和最不稳定 CAPE 的均方根差大约 10 J/kg。最佳气块和 ``unipost.exe`` 的剩余差别是预期的，也是有意保留的：upp-py 直接解湿绝热和 Bolton θe，UPP 的 ``CALCAPE`` 用查表。θe 几乎相同时两边可能选到相邻一层。大约 0.2% 的点把 RMSE 拉到约 100 J/kg，99 分位绝对差约 20 J/kg。

这三类 CIN 的平均偏差约 1–2 J/kg，最大绝对差约 680–750 J/kg。离群点上直接求解没有正浮力层，按 ``CALCAPE.f`` 把 CIN 留在 0；气块在约 350–400 hPa 只比环境冷 0.01–0.04 K。参考场在同一点有约 1–9 J/kg 的 CAPE 和数百 J/kg 的负 CIN，是查表把这一层报成略暖、从而打开从抬升凝结高度到假平衡高度的积分。浮力符号以直接求解为准，求解器不改为对齐这几百 J/kg。例子：(173, 51) 42.688°N、124.422°W，地面 CIN 0 对 −748；(151, 55) 42.076°N、124.253°W，0 对 −746；(92, 53) 40.432°N、124.282°W，混合层 0 对 −685。

0–1 km 和 0–3 km 螺旋度与 ``CALHEL.f`` 一致。可降水量按 ``CALPW.f`` 只积分水汽。这个时次的累积降水两边都是 0。

合成反射率按 NCAR + Thompson 的写出路径，取模式 ``REFL_10CM`` 的柱最大。该变量在这个文件里全是 0，所以两边都是 0，相关系数没有定义。unipost 日志里的 ``dbze_max`` 约 1 dBZ，说明诊断算过，但没有写进 REFC 那条记录。1 km 反射率的 Stoelinga 公式有 0.1 mm\ :sup:`6` m\ :sup:`−3` 的对流下限，裸算约 −10 dBZ；写出的 REFD 把晴空和低于 0 dBZ 的值记为 0，这是有意的。这份 GRIB 的 1 km 场整场为 0，两边因此一致。

整场为常数时相关系数没有定义，表里是 nan。

这张表没有能见度、结冰类别、气溶胶和亮温。该时次是起报时刻。这个 wrfout 是 ``mp_physics=8``，没有 ``PREC_ACC_*`` 和 ``smoke``。``QNWFA`` 只在 ``imp_physics==28`` 时读。VIS 和 TIPD 已按 ``CALVIS`` / ``CALICING`` 实现，没有和这份 GRIB 逐场对过。亮温不在 ``WRFPRS_d01.00`` 里。这台环境没有编出 ``unipost.exe``，所以没有 ABI 通道的偏差表。经过和原因写在 ``docs/validation.md``。
