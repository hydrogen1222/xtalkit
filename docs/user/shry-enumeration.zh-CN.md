# 安全运行一次 SHRY 枚举

这是一份生产计算手册。第一次处理真实的大型占位无序体系前，请完整读一遍。

SHRY 与 enumlib 都能生成对称不等价的有序构型。中小原胞可直接使用
`xtalkit enumerate`；当你需要在生成前知道精确数量，并希望保存每一步的
检查记录时，使用 SHRY：

```text
输入 CIF
   │
   ├─ prepare  → 补全并检查部分占据 CIF
   ├─ count    → 只计算不等价构型数
   ├─ enum     → 生成原始与清理后的有序结构
   ├─ verify   → 检查数量、化学式、空位和重复结构
   └─ postprocess（可选）→ 排序或生成 tblite 作业
```

## 1. 安装

SHRY 1.1.8 依赖旧版 pymatgen，必须与 xtalkit 隔离：

```bash
uv python install 3.11
uv tool install --python 3.11 'shry==1.1.8'
shry --help
uv tool list
```

不要在本仓库执行 `uv pip install shry`。如果命令不在 `PATH` 中：

```bash
export XTALKIT_SHRY_CMD="$(command -v shry)"
```

xtalkit 一侧的独立对称性与去重检查使用当前 pymatgen/spglib：

```bash
uv sync --extra enumerate
```

除非还要做 enumlib 交叉检查，否则不必编译 enumlib。

## 2. 先确认输入代表什么

开始前必须知道：

1. CIF 是原胞还是常规胞；
2. 目标母相的空间群与 setting；
3. 哪些 CIF 行属于同一个无序晶体学位点；
4. 乘上缩放矩阵并去除空位后的准确整胞化学式；
5. 每个部分占据能否对应整数个原子。

不要从文件名猜这些信息，应对照 CIF 标签、精修表和化学计量。

### F 心和 I 心常规胞

SHRY 不会自动转原胞。常规胞与原胞的组合数可能相差许多个数量级：从 48
个位点中选 24 个约有 (10^{13}) 种原始组合，从 12 个位点中选 6 个只有
924 种。

优先从产生该结构的精修或晶体学软件导出原胞。必须自行转换时，可用下面的
pymatgen 脚本写出“列出原胞全部位点的 P1 CIF”：

```python
import gemmi
from pymatgen.core import Structure
from pymatgen.io.cif import CifWriter
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

source = Structure.from_file("conventional.cif")
analyzer = SpacegroupAnalyzer(source, symprec=0.01, angle_tolerance=5)
primitive = analyzer.get_primitive_standard_structure()

# Pymatgen 会让同一无序位点的不同元素行共用标签；写出后将每一行改成
# 唯一标签，使最终 CIF 符合 atom-site 标识要求。
CifWriter(primitive, symprec=None).write_file("primitive_p1_raw.cif")
document = gemmi.cif.read_file("primitive_p1_raw.cif")
labels = document.sole_block().find_loop("_atom_site_label")
for index in range(len(labels)):
    labels[index] = f"site{index + 1}"
document.write_file("primitive_p1.cif")

print("常规胞位点数：", len(source))
print("原胞位点数：", len(primitive))
print("原胞成分：", primitive.composition)
```

转换后必须人工核对位点数、成分、晶格和部分占据，再删除或归档中间文件
`primitive_p1_raw.cif`。下一步总是明确给出 `--parent-spacegroup`；如果
独立检测不能恢复目标空间群，xtalkit 会拒绝最终输入。

不要把常规 F/I 胞的对称操作直接套在任意原胞矢量上。列出原胞全部位点的
P1 文件虽然不简洁，但含义明确。

## 3. 准备输入

`--target-formula` 是 SHRY 实际枚举晶胞的**准确化学式**：已经乘上缩放
矩阵的行列式，并且不含空位。程序不会自动把它约成最简式。

```bash
uv run xtalkit shry prepare primitive_p1.cif \
  --out ready.cif \
  --vacancy-symbol X \
  --parent-spacegroup 216 \
  --target-formula Li6PS5Cl \
  --scaling-matrix 1 1 1
```

`prepare` 会：

- 用显式 `X` 行补足未占满的位点；
- 拒绝超范围占据率及同一位点占据率之和大于 1；
- 检查整数原子数与准确目标化学式；
- 独立检测母相对称性并记录位点轨道。

输出为：

```text
ready.cif
ready.cif.manifest.json
ready.cif.orbit_grouping.json
```

继续前检查 manifest 中的 `parent_spacegroup_detected`、
`supercell_matrix`、`site_orbits` 和
`target_formula_after_removing_vacancy`。

若精修小数确实需要理想化，必须显式记录：

```bash
uv run xtalkit shry prepare input.cif \
  --out ready.cif \
  --set-occupancy "Li1:11/16 M1:Ge1/2P1/2" \
  --parent-spacegroup 137 \
  --target-formula Li20Ge2P4S24
```

冒号前的标签必须存在于输入 CIF。程序不做静默四舍五入。

## 4. 先计数

```bash
uv run xtalkit shry count ready.cif \
  --scaling-matrix 1 1 1 \
  --symprec 0.01 --angle-tolerance 5 --atol 1e-5 \
  --out count.json

python -m json.tool count.json
```

`count_only_result` 是将要生成的不等价结构数；`total_combinations`（若 SHRY
提供）是尚未除去对称等价关系的组合数。同时核对输入哈希、缩放矩阵、容差
和后端版本。

命令结束不等于适合继续。先估算 `count_only_result` 对应的磁盘和后处理
成本。数量异常大时应回查晶胞与无序模型，而不是立即生成。

## 5. 在新目录中生成

把计数结果原样填入 `--expect-count`：

```bash
uv run xtalkit shry enum ready.cif \
  --scaling-matrix 1 1 1 \
  --expect-count <COUNT_FROM_JSON> \
  --remove-vacancy X \
  --target-formula Li6PS5Cl \
  --write-cif --write-poscar \
  --out SHRY_run
```

`SHRY_run` 必须不存在或为空。成功目录包括：

```text
SHRY_run/
├── input/       输入副本、命令、SHRY 环境和 mod-only 结果
├── raw_shry/    SHRY 原始 CIF，仍含空位伪元素
├── clean_cif/   去除空位后的有序 CIF
├── poscar/      可选 POSCAR
├── checks/      验证与排序报告
└── manifest.json
```

正式生成前，xtalkit 会执行 SHRY 的 `--mod-only` 路径，并比较其展开成分与
准备输入。严格模式下，任何化学变化都会停止生成。

失败目录会保留，便于诊断。查看 `manifest.json` 的 `failure_stage`、输入
目录中的命令和原始输出。重试时使用新目录，不要混入旧结果。

## 6. 验证

```bash
uv run xtalkit shry verify SHRY_run \
  --target-formula Li6PS5Cl
```

默认检查生成数量、准确化学式、占据率、残留空位以及经
`StructureMatcher` 确认的重复结构。

可再做容差敏感性扫描：

```bash
uv run xtalkit shry verify SHRY_run \
  --target-formula Li6PS5Cl \
  --symprec-list 1e-4 1e-3 1e-2
```

不同容差的计数必须稳定。若发生变化，应检查近似对称坐标，而不是对结果
取平均。

### 简并度限制

`--write-degeneracy` 和 `--check-degeneracy` 只适用于缩放矩阵行列式为 1，
且 CIF 声明空间群与独立检测空间群一致的情况。真正超胞还需要母胞群作用
中没有表达的平移操作，因此程序会拒绝声称严格简并度。

列出全部原胞位点的 P1 文件仍可被识别并枚举，但不满足“声明群等于检测群”
这一条件；除非有正确声明的表示，否则不要启用简并度。

## 7. 可选后处理

按最短 Li–Li 距离排序：

```bash
uv run xtalkit shry postprocess SHRY_run \
  --shortest-distance --pair Li Li
```

按电中性的有序 Ewald 能排序：

```bash
uv run xtalkit shry postprocess SHRY_run \
  --ewald --ewald-charges Li:1 P:5 S:-2 Cl:-1
```

生成周期 tblite 作业与 Slurm 数组：

```bash
uv run xtalkit shry postprocess SHRY_run \
  --write-tblite --charge 0 --uhf 0
uv run xtalkit shry postprocess SHRY_run --write-slurm
```

按主 README 安装 tblite。执行前逐项核对 `run.sh` 中的总电荷与未成对电子
数；默认值只代表中性单重态。

## 8. 论文计算应保留什么

- 原始 CIF 与准备后的 CIF；
- 目标化学式和缩放矩阵；
- 母相空间群、setting、检测到的 Hall 符号；
- SHRY 版本和隔离环境包清单；
- count JSON、成功 manifest、verify 报告与容差扫描；
- 执行时的 xtalkit 代码修订。

计数与生成数必须相等；每个清理后的结构都必须具有目标化学式、没有空位
伪元素，也没有经结构匹配确认的重复项。

## 9. 常见失败

**检测群与 `--parent-spacegroup` 不同**
: 检查晶胞、setting、坐标和对称容差。未确认物理意义前，不要直接接受
  检测到的新空间群。

**占据率无法整数化**
: 当前轨道多重度与超胞行列式不能产生整数原子数。修正无序模型，或选择有
  晶体学依据的缩放矩阵。

**`--mod-only` 改变了化学成分**
: 比较 `input/shry_ready.cif` 与 `input/shry_mod_only.cif`，查明原因前不要
  放宽严格模式。

**生成数与 count 不一致**
: 保留失败目录，比较命令、容差、输入哈希和 SHRY 版本。不要把不同运行的
  部分结果合并。

**输出目录非空**
: 换一个新目录。这项限制用于防止旧结构被误计为本次结果。
