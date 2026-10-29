气溶胶
======

按 ``INITPOST.F`` 实际去读的变量，以及 ``MDLFLD.f`` 的单位换算。

数浓度
------

``QNCLOUD``、``QNWFA``、``QNIFA`` 只在 ``imp_physics == 28``（气溶胶感知的 Thompson）时读取。``MDLFLD.f`` 把小于 ``1e-8`` 的值写成 0。单位沿 wrfout，是 kg\ :sup:`−1`。其它微物理选项请求这些短名时，属性 ``skipped`` 写明 ``imp_physics``。Carr Fire 是 ``mp_physics=8``，不在这个分支。

质量浓度和光学厚度
------------------

``smoke`` 乘 ``PMID/(Rd*T)``，从 ug kg\ :sup:`−1` 换成 ug m\ :sup:`−3`。和 ``MDLFLD.f`` 里烟的那一行相同，不用虚温。

``TAOD5502D``、``AERASY2D``、``AERSSA2D`` 原样写出，短名 AOTK、AERASY、AERSSA。wrfout 里没有就不建场。

沙尘
----

``DUST_1`` 到 ``DUST_5`` 在 v4.1 的 ``INITPOST.F`` 里包在 ``if(1==2)`` 中，``unipost`` 不会读。文件里如果有这些数组，这里仍按 ``MDLFLD.f`` 乘密度，单位 ug m\ :sup:`−3`。没有数组时 ``skipped`` 会写这段原因。

海盐、硫酸盐、有机碳和黑碳在 ``MDLFLD.f`` 有写出。WRF 的 ``INITPOST.F`` 没有对应的 ``VarName``。不做。

亮温在 :doc:`brightness`。
