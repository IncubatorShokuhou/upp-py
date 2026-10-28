"""
Carr Fire 个例
==============

设置环境变量 ``UPP_PY_WRFOUT`` 指向 wrfout 之后，本脚本会画 Shuell 海平面气压、
500 hPa 高度和温度、合成反射率以及地面 CAPE。没有这个文件时只画说明，
文档构建不会去下载 72 MB 的历史文件。
"""

import os

import matplotlib.pyplot as plt

path = os.environ.get("UPP_PY_WRFOUT", "")
if not path or not os.path.isfile(path):
    fig, ax = plt.subplots(figsize=(6.4, 3.2))
    ax.set_axis_off()
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.text(
        0.5,
        0.55,
        "未设置 UPP_PY_WRFOUT\n跳过 Carr Fire 绘图",
        ha="center",
        va="center",
        fontsize=14,
    )
    ax.text(
        0.5,
        0.25,
        "export UPP_PY_WRFOUT=/path/wrfout_d01_2018-07-23_21_00_00",
        ha="center",
        va="center",
        fontsize=9,
        family="monospace",
    )
    fig.tight_layout()
    plt.show()
else:
    from upp_py.config import PostConfig
    from upp_py.post import post_wrfout

    cfg = PostConfig(
        fields=("PRMSL", "HGT", "TMP", "REFC", "CAPE"),
        pressure_pa=(50000.0,),
    )
    ds = post_wrfout(path, cfg)
    lat = ds["latitude"].values
    lon = ds["longitude"].values

    def _panel(ax, field, title, cmap):
        im = ax.pcolormesh(lon, lat, field, shading="auto", cmap=cmap)
        ax.set_title(title)
        ax.set_aspect("equal")
        fig.colorbar(im, ax=ax, fraction=0.046)

    fig, axes = plt.subplots(2, 2, figsize=(9.2, 8.0))
    _panel(axes[0, 0], ds["prmsl"].values[0] / 100.0, "PRMSL (hPa)", "coolwarm")
    z500 = ds["gh"].values[0, 0]
    t500 = ds["t"].values[0, 0]
    _panel(axes[0, 1], z500, "500 hPa height (m)", "viridis")
    axes[0, 1].contour(lon, lat, t500, levels=8, colors="k", linewidths=0.4)
    _panel(axes[1, 0], ds["refc"].values[0], "REFC (dBZ)", "turbo")
    _panel(axes[1, 1], ds["cape"].sel(parcel="surface").values[0], "CAPE surface (J/kg)", "YlOrRd")
    for ax in axes.ravel():
        ax.set_xlabel("longitude")
        ax.set_ylabel("latitude")
    fig.tight_layout()
    plt.show()
