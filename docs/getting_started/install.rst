安装
====

需要 Python 3.10 及以上。核心依赖是 NumPy 和 xarray、netCDF4，以及 SciPy 与 PyYAML。没有 Fortran 或 C 扩展。

.. code-block:: bash

   pip install -e ".[dev]"

读参考 GRIB、以及自己写 GRIB2，用同一个可选依赖。算法仍是 NumPy，不调用 ``unipost.exe``。``eccodes`` 的 wheel 带 C 库。

.. code-block:: bash

   pip install -e ".[grib]"

亮温 BRTMP 用 https://github.com/IncubatorShokuhou/crtm-py ，发行名 crtm-py。``pip install -e ".[crtm]"`` 会从 git 安装它。系数另外用 ``crtm-py fetch`` 下载，见 :doc:`../user_guide/brightness`。

文档用 Sphinx 和 pydata-sphinx-theme，和 scikit-learn 同一套主题。图集默认不读取 72 MB 的 wrfout；Carr Fire 那张图只有设置了 ``UPP_PY_WRFOUT`` 才会真正画个例。

.. code-block:: bash

   pip install -e ".[docs]"
   cd docs && make html

生成的页面在 ``docs/_build/html/index.html``。Read the Docs 的配置是仓库根目录的 ``.readthedocs.yaml``，构建命令就是上面的 ``make html``。

命令行入口是 ``upp-py``：

.. code-block:: bash

   upp-py --help
   upp-py --write-example-config configs/example.yaml
