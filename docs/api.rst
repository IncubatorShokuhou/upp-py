API 参考
========

文档字符串是英文，按 numpydoc 写。下面的页面由 autosummary 从源码生成。

后处理入口
----------

.. autosummary::
   :toctree: _autosummary
   :nosignatures:

   upp_py.post_wrfout
   upp_py.write_netcdf
   upp_py.load_config
   upp_py.PostConfig

读取与插值
----------

.. autosummary::
   :toctree: _autosummary
   :nosignatures:

   upp_py.wrf.open_wrfout
   upp_py.wrf.WrfState
   upp_py.interp.interpolate_pressure_levels
   upp_py.interp.PressureLevel

水汽与海平面气压，对流和反射率
--------------------------------

.. autosummary::
   :toctree: _autosummary
   :nosignatures:

   upp_py.thermo.saturation_specific_humidity
   upp_py.thermo.relative_humidity
   upp_py.thermo.dewpoint
   upp_py.thermo.fpvsnew
   upp_py.thermo.bolton_theta_e
   upp_py.post.precipitable_water
   upp_py.mslp.shuell_slp
   upp_py.mslp.membrane_slp
   upp_py.cape.cape_surface
   upp_py.cape.cape_mixed_layer
   upp_py.cape.cape_best
   upp_py.cape.cape_most_unstable
   upp_py.cape.CapeResult
   upp_py.helicity.storm_relative_helicity
   upp_py.reflectivity.reflectivity
   upp_py.clouds.cloud_fraction
   upp_py.vorticity.absolute_vorticity
