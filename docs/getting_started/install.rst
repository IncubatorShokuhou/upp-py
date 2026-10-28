安装
====

需要 Python 3.10 及以上。核心依赖是 NumPy、xarray、netCDF4、SciPy 和 PyYAML。没有 Fortran 或 C 扩展。

.. code-block:: bash

   pip install -e ".[dev]"

和参考 GRIB 对比时再装读 GRIB 的库。它们只在验证脚本里用，不参与后处理本身：

.. code-block:: bash

   pip install -e ".[grib]"

文档用 Sphinx 和 pydata-sphinx-theme，和 scikit-learn 同一套主题。图集默认不读取 72 MB 的 wrfout；Carr Fire 那张图只有设置了 ``UPP_PY_WRFOUT`` 才会真正画个例。

.. code-block:: bash

   pip install -e ".[docs]"
   cd docs && make html

生成的页面在 ``docs/_build/html/index.html``。Read the Docs 的配置是仓库根目录的 ``.readthedocs.yaml``，构建命令就是上面的 ``make html``。

命令行入口是 ``upp-py``：

.. code-block:: bash

   upp-py --help
   upp-py --write-example-config configs/example.yaml
