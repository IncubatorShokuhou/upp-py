亮温
====

``CALRAD_WCLOUD_newcrtm.f`` 用 CRTM 算卫星亮温，GRIB 短名 BRTMP。GOES-16 是传感器 ``abi_g16``（``IGET`` 927–936），GOES-17 是 ``abi_g17``（937–946），天底 ``abi_gr`` 是 958–967。每一组都是 ABI 红外通道 7 到 16。

正演在 https://github.com/IncubatorShokuhou/crtm-py 。发行名 crtm-py，导入名 ``crtm_py``。还没有发到 PyPI。可选依赖写的是

.. code-block:: text

   crtm-py @ git+https://github.com/IncubatorShokuhou/crtm-py

.. code-block:: bash

   pip install -e ".[crtm]"

upp-py 把 wrfout 的自下而上气压、温度、比湿、水成物和地表交给 ``crtm_py.brightness_temperature``。层序翻转、水汽改成 g/kg、臭氧下限、``EFFR``、陆地类型和 ``GEO_ZENITH_ANGLE`` 都在 crtm-py 里，按 ``CALRAD_WCLOUD_newcrtm.f``。``abi_g16`` 的星下点是 0°N、75.2°W。散射打开，走 ADA。返回的亮温是 ``(ny, nx, 通道)``，这里转成 ``abi_channel`` 在前。

没装 ``crtm_py``，或者目录里没有该传感器的 ``SpcCoeff.bin`` 和 ODPS ``TauCoeff.bin``，这个场不写，原因在数据集属性 ``skipped``。没有系数时 crtm-py 会用合成灰体，那个结果不当成 BRTMP。

系数用 crtm-py 自己的命令下载，不放进 upp-py：

.. code-block:: bash

   crtm-py fetch --list
   crtm-py fetch --yes

``--yes`` 会下载 CRTM REL-2.4.0 的 fix 包，大约 3.4 GB。目录默认在 ``~/.cache/crtm-py/coefficients``，也可以设 ``CRTM_PY_COEFF_DIR``。

CI 会 ``pip install -e ".[crtm]"``，不下载系数，所以 CI 里的 BRTMP 测试不跑正演。

传感器和通道可以写在 YAML 里，默认是 GOES-16 的 7–16 通道：

.. code-block:: yaml

   fields:
     - BRTMP
   brtmp_satellite: abi_g16
   brtmp_channels: [7, 8, 9, 10, 11, 12, 13, 14, 15, 16]
