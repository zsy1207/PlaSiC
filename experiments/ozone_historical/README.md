# 外部臭氧历史实验

按上游 [PR #9](https://github.com/sunmoumou1/PlaSiC/pull/9) 的反馈，比较
1850–2014 年 CO₂＋臭氧、仅臭氧、仅 CO₂ 三种强迫下的全球平均近地表气温。
本目录提供数据转换、实验启动和温度计算代码；**真实 input4MIPs 历史实验尚未完成，不能据此声称气候响应已验证。**

## 1. 数据简述

原始臭氧必须是 input4MIPs / CMIP6 / CMIP / UReading /
UReading-CCMI-1-0 / atmos / mon / vmro3 / gn / v20160711 的四个文件：

- `vmro3_input4MIPs_ozone_CMIP_UReading-CCMI-1-0_gn_185001-189912.nc`
- `vmro3_input4MIPs_ozone_CMIP_UReading-CCMI-1-0_gn_190001-194912.nc`
- `vmro3_input4MIPs_ozone_CMIP_UReading-CCMI-1-0_gn_195001-199912.nc`
- `vmro3_input4MIPs_ozone_CMIP_UReading-CCMI-1-0_gn_200001-201412.nc`

上游指定的检索入口是 <https://esgf-node.ornl.gov/search>。Google 镜像必须核对同一
文件身份及发布方校验和，不能用 CMIP6 模型输出的臭氧替代 input4MIPs 强迫。
公开元数据记录显示变量为 `vmro3(time,plev,lat,lon)`，66 个压力层、96×144 网格，
压力单位 hPa，臭氧单位 `mole mole-1`，时间单位 `months since 1850-01-01 00:00`。
转换器记录实际输入文件的 SHA256；这不能替代与发布方校验和的比对。

CO₂ 使用仓库中的 UoM-CMIP-1-2-0 月度数据。初值取
`experiments/historical_3d_500y/initial_conditions/icA.restart`，三组使用相同初值。
配置与 `historical_3d_500y/summary.json` 一致：T21、10 层大气、16 层三维海洋、2 个 MPI rank。
本轮默认只用 A 初值，不把单成员差异解释成显著的强迫响应。

## 2. 路径结构和运行方法

依赖：PlaSiC 的 C/MPI 依赖，以及 Python 3.11+、NumPy、SciPy、netCDF4、Matplotlib。
云环境每个 shell 先执行 `source /workspace/plasic-deps/activate.sh`。
所有命令从仓库根目录运行，原始数据与输出放在仓库外。

```bash
make -C src -j4 all MPI=1 NPRO=2 NLAT=32 NLEV=10
mkdir -p /workspace/plasic-runs/reference/monthly
mpirun --bind-to none -np 2 build/T21_L10_MPI_P2_O16/plasic.x \
  --restart experiments/historical_3d_500y/initial_conditions/icA.restart \
  --reset-nstep 21622368 --steps 11680 --co2-ppm 284.316853 \
  --monthly-output /workspace/plasic-runs/reference/monthly \
  --output /workspace/plasic-runs/reference/final.restart
```

上述完整 1850 年模拟只用于生成参考地面气压和检验环境，使用原有理想化臭氧；
它不是三组真实臭氧历史实验之一。后续三组仍从同一原始 `icA.restart` 启动。

取得并核验四个原始文件后：

```bash
python experiments/ozone_historical/prepare_ozone.py \
  /workspace/plasic-data/vmro3_input4MIPs_ozone_CMIP_UReading-CCMI-1-0_gn_*.nc \
  --reference-ps /workspace/plasic-runs/reference/monthly/ps.nc \
  --output /workspace/plasic-data/ozone_historical.nc
python experiments/ozone_historical/prepare_ozone.py \
  /workspace/plasic-data/vmro3_input4MIPs_ozone_CMIP_UReading-CCMI-1-0_gn_*.nc \
  --reference-ps /workspace/plasic-runs/reference/monthly/ps.nc \
  --fixed-1850 --output /workspace/plasic-data/ozone_1850_cycle.nc
python experiments/ozone_historical/run_experiments.py \
  --historical-ozone /workspace/plasic-data/ozone_historical.nc \
  --fixed-ozone /workspace/plasic-data/ozone_1850_cycle.nc
python experiments/ozone_historical/fig1/compute_temperature.py
python experiments/ozone_historical/fig1/plot_temperature.py --final
```

`PLASIC_OZONE_RUNS` 可覆盖默认结果目录 `/workspace/plasic-runs/historical`。
实验串行启动，各占两个 MPI rank；现有输出目录会报错，避免把残留结果视为成功。
`--end-year 1850` 可用于一年试运行，正式绘图计算会拒绝不完整的 165 年序列。
每组 `run.json` 保存精确命令、输入 SHA256、时间及退出状态。
`data_output/fig1/annual_temperature.csv` 保存年平均和异常；绘图再读取 CSV。
正式图输出 `fig1_surface_temperature.png` 和可编辑文字的 PDF。

## 3. 实现算法与默认假设

- `prepare_ozone.py` 逐月读取、校验连续 1980 个月以及变量、单位、缺测值。
  经纬度双线性插值，跨经度接缝周期延拓；垂直对 log(p) 线性插值。
  源压力或纬度覆盖以外取最近端点，不做线性外推。
- 固定月度参考地面气压将 sigma 映射为 `p = sigma × ps_reference(month,lat,lon)`。
  **此方法生成规定的 sigma 臭氧场，不跟随瞬时地面气压变化，也不保证臭氧柱含量守恒。**
  这是当前原生 sigma 输入接口下可复现的默认近似；正式科学结论前应检查柱含量和
  垂直分布对重网格及参考气压的敏感性。
- CO₂＋臭氧：两种强迫均随历史变化；仅臭氧：CO₂ 固定为 284.316853 ppmv；
  仅 CO₂：重复 1850 年外部臭氧季节循环，而不是改用理想化臭氧。
- 模型沿用原来的月中时间插值与 Gregorian 日历，每天 32 步；
  1850–2014 每组共 1,928,480 步。
- 温度变量为 `tas`（近地表空气温度，不是 `ts` 地表皮肤温度）。
  用 Gaussian 求积权重做全球平均、月时间界限做年平均。
  与参考图一致，每条曲线减去自身 1850–1900 年平均；不做平滑，不制造不确定性带。

## 4. 验证、问题与风险

2026-11-23 的云环境验证：

- 串行及两 rank MPI 的 `-Werror` 构建通过。
- 外部臭氧串行/MPI 联合测试 14 项通过，无跳过。
- 数据转换测试 5 项通过，包括四个合成源文件的完整 1980 月转换、固定季节循环、
  以及把两种转换结果送入真实 PlaSiC 可执行文件。
- 温度分析/实验参数测试 2 项通过，检查 Gaussian 积分、月时长加权、残缺月份拒绝和 CO₂ 对照关系。
- 四种串行网格的冷启动、重启及不匹配网格拒绝检查通过；文档严格构建通过。
- 参考模拟完整积分 11,680 步，12 个月输出，年平均全球 `tas` 为 287.000817 K，
  当前机器耗时约 246 秒。粗略外推三组完整实验约 34 小时，实际耗时需实测。

**图 1 状态：BLOCKED。** 当前环境的网络策略拒绝 Google Storage 和 ESGF，
尚未下载并验证四份真实数据；三组 1850–2014 实验及结果图未执行。
网络域名变更已保存为云环境配置草稿，保存草稿并不会自动更新当前网络。
没有以合成数据替代真实历史数据、没有将参考年模拟冒充历史验证。
正式图要求 Arial 或 Helvetica；当前机器尚未安装，字体下载也未成功。
原始数据转换、长期稳定性、臭氧柱含量误差以及最终图的视觉审核仍待完成。

## 5. 图片简述

计划图 1：三种强迫组合下 1850–2014 年全球平均近地表气温异常，
基准期为各实验自身的 1850–1900 年。该图用于比较时间演变，单成员结果不能
区分全部内部变率，也不能单凭本图证明各强迫效应线性可加。

## 6. 版本更新

- 2026-11-23：基于上游 `f29ee8e` 整理外部臭氧实现，增加 input4MIPs 转换、
  三组实验及温度分析流程；记录已执行测试和真实数据下载阻塞。
