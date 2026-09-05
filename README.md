# GEOHg

> 研究生期间 **"全球野火（开放生物质燃烧）汞排放与传输"** 研究的核心代码存档。
> 项目基于 MATLAB 实现，用于全球开放生物质燃烧（Open Biomass Burning）大气汞（Hg）排放清单的估算、不确定性量化、GEOS-Chem 模拟结果分析以及观测验证。

- **研究时段**：2010–2019（新2024）；不确定性与模拟情景分析另覆盖 2010–2019 / 2015–2024 等子时段

- **空间分辨率**：0.25°（排放清单）/ 2°*2.5°（GEOS-Chem 模拟）

- **语言**：MATLAB

***

## 模型与方法简介

排放清单采用 **质量平衡（Mass-Balance）方法**估算开放生物质燃烧的大气汞排放，基本公式如下：

```
Emission = Burned_Area × Fuel_Load × Combustion_Factor × Hg_Concentration × (1 - Ash)
```

对燃烧生物质按组分分别计算并求和，包括：

| 组分 Component                            | 说明                                                  |
| --------------------------------------- | --------------------------------------------------- |
| 叶（Leaf）、枝（Branch）、干（Bolewood）、皮（Bark） | 由地上生物量（AGB）拆分得到，受树种、树冠结构影响                      |
| 枯落物 Litterfall                          | 通过枯落物汞通量估算                                                        |
| 泥炭 Peat                                  | 泥炭地燃烧排放，与泥炭覆盖度、燃烧深度、汞含量相关                             |
| 土壤遗留汞 Soil                             | 延伸清单（`mercury_emission_new.m`）中新增，基于土壤汞排放数据并投影到燃烧网格 |

关键输入要素（对应代码内局部变量与计算主链）：

- **燃烧面积 BA**：GFED5（2002–2024，0.25°）

- **生物量燃料载荷 FL**：2010 年基准 AGB，叶/枝/干/皮由林分参数拆分，依据 NDVI/LAI 外推

- **燃烧系数 CF**：基于树冠覆盖度（MODIS VCF / Tc）与植被状况指数（VCI，由 NDVI 计算）分区估算，区分森林/灌丛/草地/农田

- **汞含量与灰分**：叶汞浓度采用机器学习全球制图产品（含观测统计数据填补），配灰分系数

- **泥炭**：覆盖度、燃烧深度（由土壤湿度驱动）与汞含量

- **GEOS-Chem 情景模拟**：对比GEOS-Chem默认GFED4s、无火灾、本研究排放清单、无北方地区火灾等情景，量化野火汞排放对大气浓度与沉降的贡献

***

## 文件说明

| 文件                           | 作用                                                                                    | 主要时段           |
| ---------------------------- | ------------------------------------------------------------------------------------- | -------------- |
| `mercury_emission.m`         | 开放生物质燃烧汞排放估算模型**核心代码**，质量平衡法逐组分计算月/年排放，输出全球、大洲、GFED 分区、土地利用、季节、典型火灾事件等统计量 | 2010–2019      |
| `Monte_Carlo.m`              | **不确定性分析**代码，对 BA、AGB、CF、汞浓度、枯落物、泥炭等参数进行蒙特卡洛模拟（10000 次），输出各栅格 P5 / P50 / P95 排放          | 2010–2019      |
| `mercury_emission_new.m`     | 在原排放计算基础上**修改数据源与处理方法**，将排放清单**延伸至 2024 年**并新增土壤（遗留）汞排放                                     | 2010–2024      |
| `GEOSChem_Output.m`          | 处理分析 GEOS-Chem 模拟的大气汞浓度与干/湿沉降，计算野火贡献、多年平均、异常与趋势                                                  | 2015–2024      |
| `Observational_validation.m` | 处理北美大气汞监测站点数据（含凋落物沉降序列），与模型模拟逐站点匹配验证                                                            | 2019–2024      |

> 说明：mercury_emission.m与Monte_Carlo.m对应已发表**EST 文章**：*Global Mercury Emissions from Open Biomass Burning Estimated Using a Mass-Balance Approach*

***

## 环境要求

- MATLAB R2020b 或更高（建议使用支持 `readgeoraster` / `geotiffwrite` 的版本）

- 常用工具箱：

  - **Mapping Toolbox**：地图栅格读写与投影

  - **Statistics and Machine Learning Toolbox**：`normrnd` / `unifrnd` / `lognrnd` / `normcdf` / `norminv` / `prctile`

  - **Image Processing Toolbox**：`imresize` 等

- 原生 **netCDF** 支持：`ncinfo` / `ncread` / `netcdf`（读取 GEOS-Chem 及 SaaS 产品）

***


- 本仓库存放代码存档，不包含大型输入/输出数据（数据文件需按结构自行组织和特殊处理），无法用于实践

