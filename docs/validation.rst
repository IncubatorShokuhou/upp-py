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

2 m 的温度、露点和相对湿度，加上 10 m 风，以及 250 hPa、500 hPa 的高度和温度，差在 GRIB 打包的量级。850 hPa 有少数点在地下外推里，温度最大差约 3 K，高度约 11 m，相关系数仍高于 0.999。

地面、混合层和最不稳定 CAPE 的均方根差大约 10 J/kg。最佳气块和 ``unipost.exe`` 的剩余差别是预期的，也是有意保留的：upp-py 直接解湿绝热和 Bolton θe，UPP 的 ``CALCAPE`` 用查表。θe 几乎相同时两边可能选到相邻一层。大约 0.2% 的点把 RMSE 拉到约 100 J/kg，99 分位绝对差约 20 J/kg。

0–1 km 和 0–3 km 螺旋度与 ``CALHEL.f`` 一致。可降水量按 ``CALPW.f`` 只积分水汽。这个时次的累积降水两边都是 0。

合成反射率按 NCAR + Thompson 的写出路径，取模式 ``REFL_10CM`` 的柱最大。该变量在这个文件里全是 0，所以两边都是 0，相关系数没有定义。1 km 反射率用的是 Stoelinga 诊断，晴空约 −10 dBZ；参考 GRIB 整场为 0，bias 因此是 −10。unipost 日志里的 ``dbze_max`` 约 1 dBZ，说明诊断算过，但没有写进 REFC 那条记录。

整场为常数时相关系数没有定义，表里是 nan。
