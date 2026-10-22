# RTS-GMLC 原始数据来源

- 仓库：https://github.com/GridMod/RTS-GMLC
- 下载提交：`3ece0d3725c844056132393ee252b3083dd4eab4`（2026-10-21由GitHub API `repos/GridMod/RTS-GMLC/commits/master`读取）
- 下载地址：`https://raw.githubusercontent.com/GridMod/RTS-GMLC/<提交>/RTS_Data/<路径>`，文件未修改。
- 许可证：该提交的仓库根目录无LICENSE文件，使用条件以仓库README及引用要求为准【文献待核查】。引用：Barrows C. et al., The IEEE Reliability Test System: A Proposed 2019 Update, IEEE Trans. Power Syst., 2020【文献待核查】。

| 文件 | 仓库路径 | SHA256 |
| --- | --- | --- |
| DAY_AHEAD_wind.csv | RTS_Data/timeseries_data_files/WIND/DAY_AHEAD_wind.csv | 6a1a8dc7d10a518523b3b319902ecc1ca1c400223832b26c6edb3a6e69d01dbc |
| REAL_TIME_wind.csv | RTS_Data/timeseries_data_files/WIND/REAL_TIME_wind.csv | fc88ecb030cec9e60a796dae9c6139855640c0fca1aa7cc76e821de24d38efb4 |
| DAY_AHEAD_pv.csv | RTS_Data/timeseries_data_files/PV/DAY_AHEAD_pv.csv | bfede6e558df5ea0f244b6326940a4ee0b95138643aa8a062897c67134c9c185 |
| REAL_TIME_pv.csv | RTS_Data/timeseries_data_files/PV/REAL_TIME_pv.csv | ae45e0114d88302520acae1c7bb278e3fdfd2691a932c13667f60a27145e0836 |
| DAY_AHEAD_regional_Load.csv | RTS_Data/timeseries_data_files/Load/DAY_AHEAD_regional_Load.csv | 7a9470d32d49068a91334cb36db54cceb2feb5cb1f702b0fa0847af8bac6cf21 |
| REAL_TIME_regional_Load.csv | RTS_Data/timeseries_data_files/Load/REAL_TIME_regional_Load.csv | 462888ae357239ec63bd22b96863900259a9849743df0ca6711ceb5a6d53da57 |
| gen.csv | RTS_Data/SourceData/gen.csv | 988466f29132b73739de60c9204dd4a2a9ceb0adf572e5966c086611272f4068 |
| bus.csv | RTS_Data/SourceData/bus.csv | cec3f776222812d43eaeaf7ac85b577719dc565f354a04e120c12a647ddb710e |

## 单位与时间约定

- 年份2020（闰年）。列`Year, Month, Day, Period`给出时段。
- 日前文件：逐小时，每日`Period`为1–24，共8784个时段；时段`h`表示当日`h-1`时至`h`时。
- 实时文件：5分钟，每日`Period`为1–288，共105408个时段；第`h`小时由`Period` `12(h-1)+1`至`12h`的平均值得到。
- 风电、光伏单位为MW（场站有功出力）；区域负荷单位为MW，列`1, 2, 3`为RTS区域1、2、3的总负荷。
- 风电4个场站`309_WIND_1, 317_WIND_1, 303_WIND_1, 122_WIND_1`，装机合计2507.9 MW；集中式光伏25座（`gen.csv`中`Unit Type = PV`），装机合计1554.5 MW。装机取`gen.csv`的`PMax MW`。
- 实时区域负荷与日前区域负荷之比在全年为常数（区域1、2、3分别约0.920、0.998、0.994，逐时标准差约2×10⁻⁵），该数据不含随机负荷预测误差。
