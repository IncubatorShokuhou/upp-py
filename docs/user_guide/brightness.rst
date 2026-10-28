亮温
====

``CALRAD_WCLOUD_newcrtm.f`` 用 CRTM 算卫星亮温，GRIB 短名 BRTMP。GOES-16 是传感器 ``abi_g16``（``IGET`` 927–936），GOES-17 是 ``abi_g17``（937–946），天底 ``abi_gr`` 是 958–967。每一组都是 ABI 红外通道 7 到 16。

正演和系数不在本仓库。独立包名 ``upp-crtm``，导入名 ``upp_crtm``。装上之后，场列表里写 ``BRTMP`` 就会调用 ``upp_crtm.brightness_temperature``。函数要接收 ``state``、``satellite``、``channels``，返回形状 ``(通道数, ny, nx)`` 的开尔文数组。系数或 CRTM 库缺失时，那个包应抛出 ``RuntimeError``，原因会进数据集属性 ``skipped``，不写出 ``brtmp``。

没装这个包时同样只记 ``skipped``，后处理其余的场照常写。不用皮肤温度代替亮温。

.. code-block:: bash

   pip install -e ".[crtm]"

``upp-crtm`` 还没发布时，这条安装会失败。CI 不装这个可选依赖。包的仓库单独建，系数也放在那边。

传感器和通道可以写在 YAML 里，默认是 GOES-16 的 7–16 通道：

.. code-block:: yaml

   fields:
     - BRTMP
   brtmp_satellite: abi_g16
   brtmp_channels: [7, 8, 9, 10, 11, 12, 13, 14, 15, 16]
