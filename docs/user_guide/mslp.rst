海平面气压
==========

这份 grib2 控制文件里有两个海平面气压，不要混：

.. list-table::
   :header-rows: 1

   * - 短名
     - 可用场编号
     - 子程序
     - 做法
   * - ``PRMSL``
     - 105
     - ``NGMSLP.f``
     - Shuell 订正的静力外推
   * - ``MSLET``
     - 23
     - ``SLP_new.f`` 的 ``MEMSLP``
     - Mesinger 膜法松弛

unipost 日志里的 ``CALLING MEMSLP`` 指的是后者。WRF 默认控制文件把 Shuell 结果写成 ``PRMSL``。Carr Fire 的参考 GRIB 里两条都有，验证时是分开比的。

Shuell（PRMSL）
---------------

从最低质量层的虚温出发，用 ``\gamma = 6.5\,\mathrm{K\,km^{-1}}`` 推到地面和海平面，得到 ``\tau = T_v R_d / g``。若海平面的 ``\tau`` 超过临界值 ``\tau_{cr} = (R_d/g)\,290.66``，就把它收到临界值上；地面和海平面都超了，再按 ``0.005\,g/R_d`` 做二次订正。然后

.. math::

   p_{sl} = p_{sfc}\,\exp(z_{sfc} / \bar\tau), \quad
   \bar\tau = \tfrac12(\tau_{sfc} + \tau_{sl})

``|z_{sfc} g| < 1`` 的点直接用地面气压。2006 年的修改把地面 ``\tau`` 改成从质量层高度而不是从地形高度往下推，代码里跟着这次修改。

膜法（MSLET）
-------------

``MEMSLP`` 在标准气压层上做。某一层的气压大于地面气压，就算地下。先找还有地下点的最高一层，对地下点做 100 次 8 点 Jacobi 松弛，过度松弛系数 ``OVERRC = 1.10``。和地下点相邻的地上点固定成虚温，作为边界。松弛始终用上一轮的值，这是 UPP 为了单进程和多进程结果一致而写的。

然后从每个点最低的地上气压层往下积分静力方程。位势穿过 0 时

.. math::

   p_{sl} = p_1 \exp\big(\Phi_1 / (R_d T)\big)

若积分到 1000 hPa 位势还是正的，再用 ``-6.6\times 10^{-4}`` 的递减率从剩下的位势外推。地面气压已经高于 1000 hPa 的点不进松弛名单，直接用地面附近一层的标准递减率外推。Carr Fire 里这类低地点和参考场的差在 10 Pa 以内。

.. code-block:: python

   from upp_py.mslp import membrane_slp, shuell_slp
   from upp_py.constants import SPL_PA
   from upp_py.interp import interpolate_pressure_levels

   prmsl = shuell_slp(state)
   levels = interpolate_pressure_levels(state, list(SPL_PA))
   mslet = membrane_slp(state, levels)

和 UPP 的差别
-------------

膜法在区域边缘不更新（Fortran 的 ``I = 2 .. IM-1``）。高山上松弛只传播约一百个格点，和 UPP 的 100 次迭代一致。参考场上 MSLET 的均方根差约 6 Pa，主要是 GRIB 打包和边缘。Shuell 的 ``PRMSL`` 同样在这个量级。
