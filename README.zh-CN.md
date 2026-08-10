# xtalkit：晶体 Wyckoff 与占位无序工具箱

[English README](README.md) · [SHRY 专题手册](docs/user/shry-enumeration.zh-CN.md)

xtalkit 用来处理几类容易出错、又经常需要重复操作的晶体学任务：

- 查询空间群及 Wyckoff 位置；
- 在现有 CIF 中标出 Wyckoff 位置，或生成纯 Wyckoff 骨架；
- 根据精修表中的空间群、晶胞和原子位置生成 CIF；
- 将部分占据结构枚举成有序构型；
- 按 Ewald 静电模型对一批可比结构排序。

它不是结构精修程序，也不会替你决定氧化态、电子自旋、合理超胞或 DFT
参数。涉及这些物理选择时，程序会要求你明确给出信息，或直接拒绝继续。

## 先从这里开始

在仓库目录中安装核心功能：

```bash
git clone https://github.com/hydrogen1222/xtalkit.git
cd xtalkit
uv sync
uv run xtalkit --version
```

用仓库自带的小文件走一遍真实流程：

```bash
# 看懂 216 号空间群有哪些 Wyckoff 位置
uv run xtalkit info --sg 216

# 在测试 CIF 中标出 4a 和 4c；结果写到 quickstart_marked.cif
uv run xtalkit mark tests/fixtures/simple.cif \
  --sg 216 --wyckoff 4a,4c -o quickstart_marked
```

`quickstart_marked.cif` 可以用 VESTA 打开。它带有
`visualization_only` 标记，只用于观察位置，不能直接拿去算能量或做定量成分分析。

不想记参数时，直接运行：

```bash
uv run xtalkit
```

在交互式终端中，这会打开终端界面。重定向输出或在 CI 中运行时则显示帮助，
不会等待输入。界面覆盖空间群查询、标记、建模、骨架、enumlib
枚举和 Ewald 排序。SHRY 需要逐阶段检查计数与清单，因此只提供命令行流程。

## 我该用哪个命令？

| 你要做的事 | 命令 | 典型输入 | 主要输出 |
|---|---|---|---|
| 查询一个空间群 | `info` | 空间群号 | Wyckoff 表 |
| 在 CIF 上显示某些位置 | `mark` | CIF、空间群号、位置标签 | 可视化 CIF/XYZ |
| 只看空间群的位置骨架 | `skeleton` | 空间群号、晶胞 | 可视化 CIF/XYZ |
| 从精修参数建立物理结构 | `build` | 空间群、晶胞、原子位点 | 物理 CIF/XYZ |
| 枚举中小型占位无序 | `enumerate` | 部分占据 CIF | 有序 CIF/XYZ 集合 |
| 枚举较大的占位无序 | `shry` | 部分占据 CIF | 分阶段清单与有序结构 |
| 比较一批有序结构的静电能 | `ewald` | CIF/POSCAR 文件夹、氧化态 | 排名 CSV |
| 检查内置数据库 | `fetch` | 无 | 230 个空间群完整性结果 |

每个命令都有简短而完整的参数说明：

```bash
uv run xtalkit --help
uv run xtalkit build --help
uv run xtalkit shry enum --help
```

## 按需安装

核心安装只包含 gemmi 和 Rich，足够使用 `info`、`mark`、`skeleton`、
`build`、`fetch` 和终端界面。

### enumlib 与 Ewald

`enumerate` 和 `ewald` 需要当前版 pymatgen/spglib；`enumerate` 还需要
Fortran enumlib：

```bash
uv sync --extra enumerate
bash scripts/build_enumlib.sh
```

编译前需要 `gfortran`、`make` 和 `git`。Debian/Ubuntu 可运行：

```bash
sudo apt install gfortran make git
```

构建脚本默认单进程，避免在共享计算节点上抢占 CPU。它把二进制装到
`~/.local/share/xtalkit/bin/`，xtalkit 会自动发现，无需修改 `PATH`。
当前默认固定到已审计提交：

```text
67217859c1cfe6e1d8cfe25b76a1d14a28c6030e
```

检查实际安装版本并跑一个小测试：

```bash
cat ~/.local/share/xtalkit/bin/enumlib-commit.txt
uv run xtalkit enumerate tests/fixtures/disordered_binary.cif \
  --max-cell-size 2 --output-dir enum_smoke
# 应得到 3 个 Au/Cu 有序 CIF
```

只有在你明确要换源码或增加编译并发时，才设置
`ENUMLIB_REF=<提交或标签>`、`XTALKIT_BUILD_JOBS=<正整数>`。

### SHRY

SHRY 1.1.8 依赖旧版 pymatgen，不能装进 xtalkit 的环境。用 Python 3.11
单独安装：

```bash
uv python install 3.11
uv tool install --python 3.11 'shry==1.1.8'
shry --help
uv tool list
```

通常 `shry` 会直接进入 `PATH`。若集群环境找不到它，再设置：

```bash
export XTALKIT_SHRY_CMD="$(command -v shry)"
```

不要执行 `uv pip install shry`：它会与 xtalkit 的新 pymatgen 发生版本冲突。
严格的对称性和去重检查仍需在 xtalkit 环境运行一次
`uv sync --extra enumerate`，但不要求编译 enumlib。

### tblite（可选）

xtalkit 只生成 tblite 作业文件；真正执行周期 xTB 单点需要外部
[tblite](https://tblite.readthedocs.io/en/stable/installation.html)。推荐放在
独立的 conda-forge 前缀：

```bash
# 尚无 micromamba 时，先按官方说明安装：
# https://mamba.readthedocs.io/en/stable/installation/micromamba-installation.html
MAMBA_ROOT_PREFIX="$HOME/.local/share/micromamba" \
  micromamba create --yes \
  --prefix "$HOME/.local/share/xtalkit/tblite-0.7.0" \
  --channel conda-forge --strict-channel-priority 'tblite=0.7.0'
mkdir -p "$HOME/.local/bin"
ln -sfn "$HOME/.local/share/xtalkit/tblite-0.7.0/bin/tblite" \
  "$HOME/.local/bin/tblite"
tblite --version
```

确保 `~/.local/bin` 在 `PATH` 中。执行作业前必须按实际体系检查总电荷和
未成对电子数，不能默认认为 `--charge 0 --spin 0` 总是正确。

## 常用任务

### 1. 查询空间群

```bash
uv run xtalkit info --sg 225
```

输出包括 Hermann–Mauguin 符号、晶系、默认示意晶胞以及每个 Wyckoff
位置的多重度、位点对称性和代表坐标。默认晶胞只用于画骨架，不是任何
真实材料的实验晶格常数。

### 2. 在现有 CIF 中标记 Wyckoff 位置

```bash
uv run xtalkit mark input.cif \
  --sg 216 \
  --wyckoff 4a,16e,48h \
  --mode overlay \
  --format cif \
  -o input_marked
```

`overlay` 保留原子并加入虚拟原子；`replace` 只替换匹配位置上的元素与
标签，保留原坐标和占据率。两种输出都是可视化文件，不是新的物理模型。

常用选项：

- `--wyckoff all`：标出该空间群的全部位置；
- `--tol 0.01`：Wyckoff 约束残差容限，分数坐标单位；
- `--map 4a:Xe,16e:Kr`：指定用于显示的虚拟元素；
- `--format cif,xyz`：同时写 CIF 和 XYZ；
- `--offset`：纯显示偏移。非零时原子会离开精确 Wyckoff 位置。

程序会检查 `--sg` 是否与 CIF 声明的空间群一致。若不一致，请先查清 CIF
采用的 setting，而不是为了让命令通过而随意改空间群号。

### 3. 生成 Wyckoff 骨架

```bash
uv run xtalkit skeleton \
  --sg 216 --wyckoff all \
  --cell "9.86 9.86 9.86 90 90 90" \
  -o SG216_skeleton
```

骨架只含代表不同位置的虚拟原子，用于教学、选位和 VESTA 对照。若省略
`--cell`，程序使用符合晶系的示意晶胞；不要把它解释为材料参数。

### 4. 从精修参数建立 CIF

如果精修表给出 Wyckoff 标签和自由参数，用 `--atom`。下面生成 NaCl：

```bash
uv run xtalkit build \
  --sg 225 \
  --cell "5.64 5.64 5.64 90 90 90" \
  --atom "Na 4a" \
  --atom "Cl 4b" \
  -o NaCl
```

例如 `16e` 的坐标模板是 `(x,x,x)`，则 `"Li 16e 0.30"` 中的 `0.30`
是自由参数；若还要给占据率，写成 `"Li 16e 0.30 0.6875"`。

如果论文或精修表直接给分数坐标，优先使用 `--atom-frac`：

```bash
uv run xtalkit build \
  --sg 137 \
  --cell "8.69407 8.69407 12.5994 90 90 90" \
  --atom-frac "Li 0.2563 0.2718 0.1832 0.6875" \
  --atom-frac "P 0 0 0.5 1" \
  -o refinement_model
```

程序会报告识别到的 Wyckoff 位置，便于发现坐标、setting 或空间群抄错。
同一位置的不同元素可用多条原子记录表示部分占据；它们的占据率之和不能
超过 1。占据率为 0 的 `--atom-frac` 行会明确跳过。

重复建模时可使用 JSON：

```json
{
  "sg": 225,
  "cell": {
    "a": 5.64, "b": 5.64, "c": 5.64,
    "alpha": 90, "beta": 90, "gamma": 90
  },
  "atoms": [
    {"element": "Na", "wyckoff": "4a"},
    {"element": "Cl", "wyckoff": "4b"}
  ]
}
```

```bash
uv run xtalkit build --spec nacl.json -o NaCl
```

`--spec` 与 `--sg`、`--cell`、`--atom`、`--atom-frac` 互斥。

### 5. 用 enumlib 枚举占位无序

先用仓库小体系确认安装：

```bash
uv run xtalkit enumerate tests/fixtures/disordered_binary.cif \
  --min-cell-size 1 --max-cell-size 2 \
  --output-dir AuCu_enum
```

重要概念：

- `min/max-cell-size` 是相对原胞的体积倍数范围，不是沿三个轴分别放大；
- xtalkit 会先转原胞，再调用 enumlib；
- `--max-structures` 只限制 enumlib 搜索完成后写多少文件，不能缩短
  `enum.x` 的搜索；
- `--jobs` 只并行生成和解析结果，`enum.x` 本身仍是单线程；
- 默认 `jobs=1`、分批写盘，内存占用最稳妥。

对新体系建议先用最小范围：

```bash
uv run xtalkit enumerate disordered.cif \
  --min-cell-size 1 --max-cell-size 1 \
  --timeout 10 --max-structures 20 \
  --output-dir trial_enum
```

若占据率在给定超胞范围内无法变成整数原子数，程序会在启动 enumlib 前
拒绝运行并给出最接近的可行占据率。不要轻易使用 `--skip-preflight`；错误
输入可能让 enumlib 消耗极大内存。

### 6. 用 SHRY 枚举较大体系

SHRY 工作流刻意分成五步：

```text
prepare → count → enum → verify → postprocess（可选）
```

先准备输入和检查成分：

```bash
uv run xtalkit shry prepare input.cif \
  --out ready.cif \
  --parent-spacegroup 216 \
  --target-formula Li6PS5Cl \
  --scaling-matrix 1 1 1
```

只计数，不生成：

```bash
uv run xtalkit shry count ready.cif \
  --scaling-matrix 1 1 1 --out count.json
```

确认 `count.json` 中的 `count_only_result` 合理后，再把它填入生成命令：

```bash
uv run xtalkit shry enum ready.cif \
  --scaling-matrix 1 1 1 \
  --expect-count <COUNT> \
  --target-formula Li6PS5Cl \
  --write-cif --write-poscar \
  --out SHRY_run

uv run xtalkit shry verify SHRY_run \
  --target-formula Li6PS5Cl
```

F 心或 I 心的大常规晶胞可能产生天文数量的组合。必须理解原胞、对称操作
和目标化学式后再运行。完整说明、输出目录和失败恢复方法见
[SHRY 专题手册](docs/user/shry-enumeration.zh-CN.md)。

### 7. Ewald 排序

Ewald 能只适用于已经有明确氧化态、完全有序且电中性的周期结构。下面按
总静电能从低到高排列同成分、同原子数的一批结构：

```bash
uv run xtalkit ewald ordered_cifs/ \
  --charges Li:1 P:5 S:-2 Cl:-1 \
  --sort asc \
  --out ewald_results/ranking.csv
```

筛选前 20 个并复制到分组目录：

```bash
uv run xtalkit ewald ordered_cifs/ \
  --charges Li:1 P:5 S:-2 Cl:-1 \
  --top-n 20 --group --group-dir ewald_results
```

默认是复制。只有明确使用 `--move` 时才移动原文件。`--guess` 可让 pymatgen
猜氧化态，但论文计算应检查猜测结果；不要把 Ewald 模型能当成弛豫后的
DFT 总能。不同成分不能放在同一次排名中；只有同成分的可公度晶胞才适合
用 `--per-atom` 比较。

## 输出文件怎么读

- `mark`、`skeleton`：带 `visualization_only` 审计标记的可视化文件；
- `build`：含空间群操作、占据率和晶胞化学式的物理 CIF；
- `enumerate`：默认写入 `<输入名>_enum/`；每个文件是一种有序构型；
- `shry enum`：`input/` 保存输入与环境，`raw_shry/` 保存原始结果，
  `clean_cif/` 保存去除空位后的 CIF，`manifest.json` 记录本次运行；
- `ewald`：完整排名写入 CSV；使用 `--group` 时另建 `selected/` 与 `rest/`。

默认输出可能覆盖同名的单个 CIF/XYZ。重要任务请总是显式给 `-o` 或新的
输出目录。SHRY 枚举目录必须是新目录或空目录，防止不同批次混在一起。

## 常见问题

### `VIRTUAL_ENV=... does not match .../.venv`

当前 shell 激活了另一个虚拟环境。通常先运行 `deactivate`，再使用
`uv run ...`。如果你确实要使用已激活环境，需明确使用 uv 的 `--active`，
不要忽略环境混用带来的依赖差异。

### 找不到 `enum.x` 或 `makestr.x`

重新运行 `bash scripts/build_enumlib.sh`。无需把安装目录加入 `PATH`；xtalkit
会从 `~/.local/share/xtalkit/bin/` 自动加载。可检查：

```bash
ls -l ~/.local/share/xtalkit/bin/enum.x \
      ~/.local/share/xtalkit/bin/makestr.x
```

### `non-integer stoichiometry`

部分占据乘以该 Wyckoff 轨道的位点数和超胞倍数后，不是整数。扩大
`--max-cell-size`，或根据原始精修与化学计量修正占据率。不要只为绕过检查
而使用 `--skip-preflight`。

### 空间群不一致

确认 CIF 中的 Hermann–Mauguin 符号、International Tables 编号、setting
和对称操作是否彼此一致。程序拒绝猜测，因为错误空间群会改变 Wyckoff
轨道、化学式和枚举数。

### SHRY 数量太大或长时间不结束

先停止生成计划，回到 `shry count`。检查是否误用了 F/I 心常规胞、是否
忘记设置正确缩放矩阵，以及部分占据是否对应目标化学式。不要直接删除
失败目录；先查看 `manifest.json` 和 `input/shry_mod_only.cif`。

### Ewald 报不电中性

检查每种元素的氧化态以及结构化学式。程序不会添加均匀补偿背景，也不会
替你修改电荷以得到一个数值。

### TUI 输错了怎么办

关键数值会要求重新输入，不会静默改成危险默认值。当前工作流若无需继续，
可按 `Ctrl+C` 退出；已经写出的文件不会自动删除。

## 开发与验证

```bash
uv sync --extra enumerate
uv run pytest
uvx ruff check xtalkit tests
uv build
```

项目要求所有 230 个空间群均可读取；涉及枚举、对称性、成分、空位、去重
和 Ewald 中性的行为都有回归测试。真实 SHRY、enumlib 或 tblite 可执行文件
存在时，相应集成测试会自动运行；不存在时会明确显示为 skipped。
