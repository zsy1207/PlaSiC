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

四个文件合计 5,927,876,651 字节。ESGF 目录中 ORNL 和 Argonne 两个副本给出的 SHA256 一致，
已记录于 [SHA256SUMS](SHA256SUMS)。从仓库根目录执行以下命令即可从 ORNL 下载并核验。

```bash
manifest="$PWD/experiments/ozone_historical/SHA256SUMS"
base=https://esgf-node.ornl.gov/thredds/fileServer/user_pub_work/input4MIPs/CMIP6/CMIP/UReading/UReading-CCMI-1-0/atmos/mon/vmro3/gn/v20160711
mkdir -p /workspace/plasic-data
(cd /workspace/plasic-data && while read -r sha name; do
  curl --fail --location --retry 3 --continue-at - --output "$name" "$base/$name" || exit 1
done < "$manifest" && sha256sum --check "$manifest")
```

变量 `vmro3(time,plev,lat,lon)` 存储 96×144 网格、66 个压力层的臭氧摩尔分数。
压力单位为 hPa，时间轴为相对于 1850 年 1 月的月数。转换程序再次核验输入 SHA256，
并在输出中记录源数据和参考气压的校验和。原始数据由 University of Reading 的 CCMI 团队提供。
[数据集永久标识](https://hdl.handle.net/21.14100/ab27892c-eaa5-3430-a3c5-f9b0549a556f) · [ESGF 检索入口](https://esgf-node.ornl.gov/search)

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
每组使用两个 MPI rank。可用 `--jobs 2` 在有四个 CPU 的机器上同时运行两组。
每个日历年结束时先核对 restart 时间、气温月界及强迫来源，再保存为 `<case>/<year>/`。
中断后以相同输入再次运行，程序核对每个已完成年份的 restart 和月度气温文件的 SHA256，再从年末断点继续。任一摘要缺失或不符都会报错。
`--end-year 1850` 可用于一年试运行，后续可继续至 2014。每段的 `run.json` 记录命令、
时间、输入、气温与重启文件的 SHA256。未通过检查的年段保留在临时目录中供诊断。
`data_output/fig1/annual_temperature.csv` 保存年平均和异常；绘图再读取 CSV。
正式图输出 `fig1_surface_temperature.png` 和可编辑文字的 PDF。

## 3. 实现算法与默认假设

- `prepare_ozone.py` 逐月读取、校验连续 1980 个月以及变量、单位、缺测值。
  水平采用周期经度和纬度双线性插值。垂直把各源层视为在其 `bounds_plev` 内均匀，
  以源层和目标层的重叠压力厚度加权求平均。在源网格下边界至局地地面间延用最底层值。
- 采用同一个参考年各月的地面气压，将每个 sigma 界面映射为 `p = sigma × ps_reference`。
  垂直重网格在这个参考气压下逐柱守恒。辐射计算使用规定的 sigma 臭氧场和模型瞬时气压，
  两种气压的差异及十层网格的垂直分辨率是解释结果时的限制。
- CO₂＋臭氧：两种强迫均随历史变化；仅臭氧：CO₂ 固定为 284.316853 ppmv；
  仅 CO₂：重复 1850 年外部臭氧季节循环，而不是改用理想化臭氧。
- 模型沿用原来的月中时间插值与 Gregorian 日历，每天 32 步；
  1850–2014 每组共 1,928,480 步。两种臭氧输入采用相同端点约定。
  1850 年 1 月上半月保持首月值，2014 年 12 月下半月保持末月值。
- 温度变量为 `tas`（近地表空气温度，不是 `ts` 地表皮肤温度）。
  用 Gaussian 求积权重做全球平均、月时间界限做年平均。
  与参考图一致，每条曲线减去自身 1850–1900 年平均；不做平滑，不制造不确定性带。

## 4. 验证、问题与风险

2026-11-23 的云环境验证：

- 串行及两 rank MPI 的 `-Werror` 构建通过。
- 外部臭氧串行/MPI 联合测试 14 项通过，无跳过。
- 数据转换测试 5 项通过，包括四个合成源文件的完整 1980 月转换、固定季节循环、
  以及把两种转换结果送入真实 PlaSiC 可执行文件。
- 温度分析和实验参数测试覆盖 Gaussian 积分、月时长加权、残缺月份拒绝、CO₂ 对照关系与断点摘要完整性。
- 四种串行网格的冷启动、重启及不匹配网格拒绝检查通过，文档严格构建通过。
  臭氧测试入口为 `make -C src test-ozone`、`make -C src test-ozone-mpi` 和 `make -C src test-data`。
  上游 `f29ee8e` 的通用 `test-runtime` 目标引用了其提交中缺少的脚本，本次四网格回归用独立验证脚本执行。
- 参考模拟完整积分 11,680 步，12 个月输出，年平均全球 `tas` 为 287.000817 K，
  当前机器耗时约 246 秒。粗略外推三组完整实验约 34 小时，实际耗时需实测。

用正式强迫资料提前试跑了七个日历窗口。1899/1900、1949/1950、1999/2000 三个资料交界各覆盖前后两个完整月。1900 年 2 月为 28 天，2000 年为 29 天，2014 年 12 月完整输出。继续进入 2015 年的第一个时间步，模型按预期报告臭氧强迫越界。

审查员独立复算了真实 2000 年四个月的臭氧柱，并对跨年连续积分和分段积分进行了对比。
方法、数值和适用范围见 [独立审查记录](REVIEW.md)。正式图需三组各 1,980 个完整月份。

## 5. 图片简述

计划图 1：三种强迫组合下 1850–2014 年全球平均近地表气温异常，
基准期为各实验自身的 1850–1900 年。该图用于比较时间演变，单成员结果不能
区分全部内部变率，也不能单凭本图证明各强迫效应线性可加。

## 6. 版本更新

- 2026-11-23：基于上游 `f29ee8e` 整理外部臭氧实现，增加 input4MIPs 转换、
  三组实验及温度分析流程。根据真实源数据的臭氧柱核算，采用压力层重叠加权，并为长积分增加年度断点。
