[English](README.md) | **中文**

# xtalkit — 晶体 Wyckoff 工具箱

在晶体结构中用虚拟原子标记 Wyckoff 位置，方便在 [VESTA](https://jp-minerals.org/vesta/en/) 中直观可视化。

## 安装

要求 Python 3.10+。

```bash
git clone <repo-url> && cd xtalkit
uv sync
uv pip install -e .
```

验证安装：

```bash
xtalkit --version
# xtalkit 0.1.0
```

## 快速开始

在一个 CIF 文件中标记两个 Wyckoff 位置：

```bash
xtalkit mark Li6PS5Cl.cif --sg 216 --wyckoff 4a,24f
# 输出：Li6PS5Cl_WYCK.cif
```

用 VESTA 打开结果 —— 虚拟原子（Xe、Kr 等）会以不同颜色标记出 Wyckoff 位置。

---

# 使用说明

xtalkit 提供两种界面：

| 界面 | 启动方式 | 适用场景 |
|------|----------|----------|
| **CLI**（命令行） | `xtalkit <command> ...` | 脚本化、批处理、一行命令 |
| **TUI**（交互式） | `xtalkit`（无参数） | 交互探索、引导式工作流 |

---

## CLI 命令参考

```
xtalkit <command> [options]
```

共七个子命令：

| 命令 | 功能 |
|------|------|
| `mark` | 在 CIF 文件中标记 Wyckoff 位置 |
| `skeleton` | 生成纯 Wyckoff 骨架（不含真实原子） |
| `info` | 查询某个空间群的 Wyckoff 位置信息 |
| `fetch` | 校验空间群数据库完整性 |
| `enumerate` | 枚举对称不等价的有序构型（需 pymatgen + enumlib） |
| `shry` | 面向大规模部分占据枚举的严格 SHRY 工作流 |
| `build` | 由精修参数生成 CIF（空间群 + 晶胞 + Wyckoff 位点 + 占据率） |

---

### `mark` — 在 CIF 中标记 Wyckoff 位置

```
xtalkit mark <input.cif> --sg <N> --wyckoff <letters> [options]
```

**必填参数：**

| 参数 | 说明 |
|------|------|
| `input.cif` | CIF 文件路径（支持相对路径） |
| `--sg N` | 空间群编号，1–230 |
| `--wyckoff L` | 要标记的 Wyckoff 字母，逗号分隔（如 `4a,24f`），或 `all` 表示全部 |

**可选参数：**

| 选项 | 默认值 | 说明 |
|------|--------|------|
| `--mode overlay` | overlay | `overlay`：保留真实原子 + 添加虚拟原子。`replace`：将匹配的真实原子替换为虚拟原子 |
| `--tol 0.5` | 0.5 | 分数坐标匹配容差 |
| `--map 4a:Xe,16e:Kr` | （自动） | 按 Wyckoff 字母自定义虚拟元素分配 |
| `--format cif` | cif | 输出格式：`cif`、`xyz`，或逗号分隔如 `cif,xyz` |
| `-o base` | `{name}_WYCK` | 输出基础路径（每种格式自动加扩展名） |

**示例：**

```bash
# 在立方 F-43m 结构中标记 4a 和 24f（overlay 模式）
xtalkit mark Li6PS5Cl.cif --sg 216 --wyckoff 4a,24f

# 标记全部 8 个 Wyckoff 位置，输出所有格式
xtalkit mark structure.cif --sg 216 --wyckoff all --format cif,xyz

# Replace 模式：把 4a 位置的真实原子替换为 Xe 虚拟原子
xtalkit mark structure.cif --sg 216 --wyckoff 4a --mode replace

# 自定义元素映射 + 紧容差
xtalkit mark NaCl.cif --sg 225 --wyckoff 4a,4b --map 4a:He,4b:Ne --tol 0.01

# 显式指定输出路径
xtalkit mark input.cif --sg 216 --wyckoff 4a -o ./output/marked
# 生成：./output/marked.cif
```

**输出文件：**

| 格式 | 文件 | 说明 |
|------|------|------|
| CIF | `{name}_WYCK.cif` | 完整 CIF，末尾追加虚拟原子。可直接用 VESTA 打开。 |
| XYZ | `{name}_WYCK.xyz` | 简单笛卡尔 XYZ（丢失晶胞信息 —— VESTA 打开时会提示输入晶胞）。 |

---

### `skeleton` — 生成 Wyckoff 骨架

生成一个**只含** Wyckoff 位置虚拟原子、不含任何真实原子的结构 —— 适合作为参考模板。

```
xtalkit skeleton --sg <N> --wyckoff <letters> [options]
```

**必填参数：**

| 参数 | 说明 |
|------|------|
| `--sg N` | 空间群编号，1–230 |
| `--wyckoff L` | Wyckoff 字母，逗号分隔或 `all` |

**可选参数：**

| 选项 | 默认值 | 说明 |
|------|--------|------|
| `--cell "a b c α β γ"` | 按晶系默认 | 自定义晶胞参数 |
| `--map 4a:Xe,...` | （自动） | 元素覆盖 |
| `--format cif` | cif | `cif`、`xyz`，或 `cif,xyz` |
| `-o base` | `SG{N}_skeleton` | 输出基础路径 |

**示例：**

```bash
# F-43m 默认立方晶胞（a=b=c=5.0）的骨架
xtalkit skeleton --sg 216 --wyckoff all

# 使用 Li6PS5Cl 真实晶胞参数生成骨架
xtalkit skeleton --sg 216 --wyckoff 4a,4c,16e,24f \
    --cell "9.85 9.85 9.85 90 90 90"

# P2₁/c（单斜）自定义晶胞骨架
xtalkit skeleton --sg 14 --wyckoff 2a,4e \
    --cell "5.5 6.3 8.1 90 108.5 90" --format cif,xyz
```

**按晶系的默认晶胞参数：**

| 晶系 | a | b | c | α | β | γ |
|------|---|---|---|---|---|---|
| 三斜 | 5.0 | 6.0 | 7.0 | 80° | 90° | 100° |
| 单斜 | 5.0 | 6.0 | 7.0 | 90° | 110° | 90° |
| 正交 | 5.0 | 6.0 | 7.0 | 90° | 90° | 90° |
| 四方 | 5.0 | 5.0 | 7.0 | 90° | 90° | 90° |
| 三方 | 5.0 | 5.0 | 8.0 | 90° | 90° | 120° |
| 六方 | 5.0 | 5.0 | 8.0 | 90° | 90° | 120° |
| 立方 | 5.0 | 5.0 | 5.0 | 90° | 90° | 90° |

**精确计算时请务必用真实晶胞参数覆盖默认值。**

---

### `build` — 由精修参数生成 CIF

```
xtalkit build --sg <N> --cell "<a b c α β γ>" --atom "<spec>" [--atom ...] [选项]
```

由 XRD 精修结果拼装标准 CIF：一个空间群、一个晶胞、一组原子位点。每个位点是一个元素落在某 Wyckoff 位置（如 `16e`）上，附带自由坐标值与可选占据率。**晶系由空间群自动推导**（不单独输入），仅用来校验晶胞参数。占据率默认 1.0；支持部分占据/混合占据，产出的无序 CIF 可再交给 `xtalkit enumerate` 枚举有序构型。

支持全部 230 个空间群（内置 Wyckoff 数据集，源自国际晶体学表 via pyxtal，并已用 gemmi 校验）。

**必选参数：**

| 参数 | 说明 |
|------|------|
| `--sg N` | 空间群号，1–230 |
| `--cell "a b c α β γ"` | 晶胞参数 |
| `--atom SPEC` | 原子规格（可重复）：`元素 wyckoff [自由坐标...] [占据率]` |
| `--spec FILE` | JSON 规格文件（替代 `--sg`/`--cell`/`--atom`） |

**原子规格格式：** `元素 wyckoff [自由坐标值] [占据率]`。自由坐标按其在位点坐标模板中出现的顺序输入；解析器知道每个 Wyckoff 位点需要几个自由值。仅当多给一个数字时，末位才被当作占据率。

| 规格 | 含义 |
|------|------|
| `Na 4a` | Na 在 4a（无自由参数），占据率 1.0 |
| `Li 16e 0.3` | Li 在 16e `(x,x,x)`，x=0.3，占据率 1.0 |
| `Li 16e 0.3 0.5` | Li 在 16e，x=0.3，占据率 0.5 |
| `S 48h 0.25 0.3` | S 在 48h `(x,x,z)`，x=0.25、z=0.3，占据率 1.0 |
| `Li 4a 0.5` + `Cu 4a 0.5` | 4a 上 Li/Cu 混合占据（无序） |

**选项：** `--format cif[,xyz]`（默认 `cif`）、`-o/--output`（输出基名；默认 `SG{N}_built`）。

**分数坐标直输模式（`--atom-frac`）** —— 手里有精修表时更省事。不用写 Wyckoff 字母 + 自由参数，直接给每个原子的最终分数坐标：`元素 x y z [占据率]`。工具自动识别每个原子落在哪条 Wyckoff 轨道（并打印出来），所以非规范代表（如 SG 137 4d 写成 `(0, 1/2, z)` 而非规范的 `(1/2, 0, 1/4+z)`，两者差一个 4 重旋转）也能直接用，无需手动算偏移。占据率为 0 的原子自动跳过。

```bash
xtalkit build --sg 137 --cell "8.694 8.694 12.599 90 90 90" \
    --atom-frac "Li 0.2563 0.2718 0.1832 0.691" \
    --atom-frac "Li 0 0.5 0.9446 1" \
    --atom-frac "S 0 0.1843 0.4103 1" \
    ...
```

**示例 —— 生成 NaCl（Fm-3m，225）：**

```bash
xtalkit build --sg 225 --cell "5.64 5.64 5.64 90 90 90" \
    --atom "Na 4a" --atom "Cl 4b" -o NaCl
# [OK] SG #225 (Fm-3m), crystal system: cubic, formula: NaCl
#      Saved to: NaCl.cif
```

**示例 —— 生成无序 Li/Cu 位点，再枚举其有序构型：**

```bash
xtalkit build --sg 225 --cell "4 4 4 90 90 90" \
    --atom "Li 4a 0.5" --atom "Cu 4a 0.5" -o AuCu_disordered
xtalkit enumerate AuCu_disordered.cif --max-cell-size 2
```

CIF 中只写不对称单元代表原子 + 空间群对称操作；VESTA、gemmi、pymatgen 会自动展开为完整晶胞。生成后 xtalkit 会打印由「多重性 × 占据率」算出的化学式（如 `NaCl`、`Li6PS5Cl`）—— 与你的精修结果对照以确认输入无误。

**JSON 规格**（`--spec`）便于复现/脚本化：

```json
{"sg": 225,
 "cell": {"a": 5.64, "b": 5.64, "c": 5.64, "alpha": 90, "beta": 90, "gamma": 90},
 "atoms": [
   {"element": "Na", "wyckoff": "4a", "occ": 1.0},
   {"element": "Cl", "wyckoff": "4b", "occ": 1.0}
 ]}
```

---

### `info` — 查询空间群信息

```
xtalkit info --sg <N>
```

**示例：**

```
$ xtalkit info --sg 216

Space Group #216: F-43m
Crystal System: cubic
Default cell: a=5.0 b=5.0 c=5.0 α=90 β=90 γ=90

Wyckoff Positions (8):
  Letter   Mult   Site Sym   Coordinates
  4a       4      -43m       0,0,0
  4b       4      -43m       1/2,1/2,1/2
  4c       4      -43m       1/4,1/4,1/4
  4d       4      -43m       3/4,3/4,3/4
  16e      16     .3m        x,x,x
  24f      24     2.mm       x,0,0
  24g      24     2.mm       x,0,1/2
  48h      48     ..m        x,x,z
```

---

### `fetch` — 校验数据库

```
xtalkit fetch
```

校验全部已内置的空间群数据是否完整。输出：

```
[OK] Space group data intact (230/230 space groups supported).
```

---

### `enumerate` — 枚举有序构型

使用 [pymatgen](https://pymatgen.org/) + [enumlib](https://github.com/msg-byu/enumlib)（Hart-Forcade 算法）枚举一个无序 CIF 的所有**对称不等价**有序构型。母相空间群对称操作相关的结构会自动合并为一个 —— 你拿到的都是真正不同的排列，无需手动去重。

本工具适用于任何**占位无序**（晶体学位点上有部分占位）的材料：二元合金（Au/Cu）、固溶体、掺杂体系，以及存在可动离子无序的电池电解质（如 argyrodite Li₆PS₅Cl、LGPS 等）。该子命令为**可选功能** —— 需要 `enumerate` uv extra 与源码编译的 enumlib 二进制（见下文 [枚举功能配置](#枚举功能配置)）。

```
xtalkit enumerate <input.cif> [options]
```

**必填参数：**

| 参数 | 说明 |
|------|------|
| `input.cif` | 含部分占位/无序占位的 CIF 文件路径（如 Au0.5/Cu0.5，或 Li0.5 + 空位） |

**可选参数：**

| 选项 | 默认值 | 说明 |
|------|--------|------|
| `--min-cell-size N` | 1 | 枚举的最小超胞尺寸 |
| `--max-cell-size N` | 2 | 枚举的最大超胞尺寸（越大→构型越多、越慢，见下文） |
| `--symm-prec TOL` | 0.1 | 空间群分析的对称容差 |
| `--vacancy-symbol S` | `X` | 内部表示空位用的 DummySpecies 符号 |
| `--output-dir DIR` | `{name}_enum/` | 枚举文件的输出目录 |
| `--max-structures N` | 不限 | 输出与后续 `makestr.x` 生成上限（不限制 `enum.x` 搜索） |
| `--timeout MIN` | 无 | enumlib 子进程超时（分钟） |
| `--format F` | cif | 输出格式：`cif` 或 `xyz` |
| `--jobs N` | 1 | 结构生成阶段的并行进程数（`0` = 自动/CPU 核数） |
| `--batch-size N` | 256 | 每批结构数 —— 限制内存峰值 |
| `--scratch-dir DIR` | 系统 temp | enumlib scratch 目录（用 `/dev/shm` 走 tmpfs） |
| `--skip-preflight` | 关 | 跳过非整数计量比预检（危险——可能 OOM） |

**示例：**

```bash
# 枚举 50/50 Au/Cu 二元体系 —— 先用 1× 超胞快速看一眼
xtalkit enumerate AuCu_disordered.cif --max-cell-size 1

# 完整 1×–2× 枚举，只保留前 5 个结构
xtalkit enumerate disordered.cif --max-cell-size 2 --max-structures 5

# 枚举一个电池电解质 CIF，自定义输出目录
xtalkit enumerate Li6PS5Cl.cif --max-cell-size 1 --output-dir ./runs/exp01/
```

> **到底枚举什么？** 你不需要在命令里指定枚举哪个元素 —— 没有这个 flag。xtalkit 从 CIF 自动判断，枚举**每一个**占位 < 1（或同一坐标上有多种物种）的位点；占位为 1 的全占位点原样固定。如果 CIF 里有多个无序位点（比如一个 Li/空位位点*外加*一个 Ge/P 混合位点），它们会被**联合**枚举 —— 每个输出结构把所有无序位点同时定死。所以 `enumerate Li6PS5Cl.cif` 枚举的是 Li 位点（Li/空位）；`enumerate AuCu_disordered.cif` 枚举的是 Au/Cu 位点。
>
> **如何选择研究对象。** 枚举完全由 CIF 决定，所以你通过"在 CIF 里写什么无序"来选择研究什么 —— 而不是靠命令行参数。只想研究 Li/空位，就把其他位点保持全占（占位 1）；只想研究 Ge/P 混合，就把 Li 保持全占。反过来，一个完全有序的 CIF（所有位点占位都是 1）没有可枚举的内容，会得到 0 个结构 —— 你得先改 CIF、把要研究的无序写进去（见下文*当 enumlib 返回 0 个结构时*）。

**工作流程：**

1. 通过 `pymatgen.core.Structure.from_file` 读取 CIF（原生支持部分占位 + mmCIF）。
2. **转换为原胞**（`Structure.get_primitive_structure()`）。这一步很关键：enumlib 不会自动找原胞，而 F 心/I 心的常规晶胞候选位点数会是原胞的 2×–4× —— 足以让枚举慢得不可行，或让 enumlib 内部数组溢出。在原胞里做候选位点数最小。
3. 对任何部分占位位点追加显式空位物种（`DummySpecies`），让 enumlib 有具体物种可放，如 `Li0.5` → `Li0.5 + X0.5`；`Au0.3/Cu0.3` → `Au0.3 + Cu0.3 + X0.4`。**占位为 1 的全占位点原样保留** —— 它们是固定 spectator，只对无序位点做枚举。
4. 调用 `pymatgen.command_line.enumlib_caller.EnumlibAdaptor`，它会 shell 调用 Fortran 二进制 `enum.x` 与 `makestr.x`。
5. 把每个不同的有序构型写为 `<basename>_<NNN>.cif`（以原胞形式输出，报告为 P1）。

#### 如何选择 `--max-cell-size`

enumlib 在原胞的**超胞**里枚举有序构型，从 `--min-cell-size` 到 `--max-cell-size`。最大尺寸越大，探索的超胞越大、找到的构型越多 —— 但数量（和耗时）可能爆炸式增长。

- **先用小尺寸。** 第一次看建议 `--max-cell-size 1`。如果母相占位在原胞里已经能整数化（如 1/2、1/3、2/3），1× 超胞就能给出全部不等价构型，几秒就跑完。
- **确有需要再放大。** 2× 超胞能揭示 1× 装不下的额外构型，但也可能产出上千个结构、跑几分钟到几小时。搭配 `--max-structures` 与 `--timeout` 控制规模。
- **占位必须可整数化。** 所选超胞里每个物种的总数必须是整数。占位 1/2 至少需要 2 个位点（2 位点原胞的 1× 即可）；1/3 需要 3× 超胞；像 0.56（= 14/25）这种值需要 25× 超胞，不现实。如果你的 CIF 里有这种值，先把它舍入到附近的简单分数（见下文）。

#### 内存、磁盘与并行

大规模枚举（上千个结构）会吃大量内存和磁盘。xtalkit 采用**分批流式**处理（不会把所有结构同时留在内存），并提供三个旋钮：

- **`--max-structures N`** 限制的是 `enum.x` 结束之后的 `makestr.x` + 解析 + 写文件阶段，不限制 `enum.x` 搜索本身。它能避免只想采样时仍生成上千个 `vasp.*`/CIF 中间文件，但 `enum.x` 仍必须先完成完整对称搜索。
- **`--batch-size N`** 每批生成/解析/写入多少个结构。内存峰值随批大小增长，而非总结构数。内存不够就调小（如 `--batch-size 32`）；想少调几次 `makestr.x` 就调大。
- **`--jobs N`** 多进程并行各批，加速结构生成（makestr + 解析 + 写文件）阶段。`--jobs 0` 自动用满 CPU。**每个 worker 会加载一份 pymatgen（约 200 MB），所以内存随 `--jobs` 增长** —— 按你的内存选合适的值。
- **`--scratch-dir /dev/shm`** 把 enumlib 的 scratch 文件（`struct_enum.out`、各批 `vasp.*`）放到 tmpfs，避免磁盘 I/O。`/dev/shm` 通常只有内存的一半大小，`struct_enum.out` 装得下再用。

> **`enum.x` 是单线程的。** 枚举核心 `enum.x` 是串行 Fortran 回溯搜索 —— `--jobs` **不能**加速它。`--jobs` 只并行其后的 makestr + 解析 + 写文件阶段。若 `enum.x` 占大头（搜索空间巨大、长时间枚举），`--jobs` 对墙钟时间提升有限；流式与 `--max-structures` 只在 `enum.x` 产出 `struct_enum.out` 之后生效。要并行 `enum.x` 本身需按晶胞大小拆分，暂未实现。

#### CIF 的分数占位是怎么工作的

在 CIF 中，同一个晶体学坐标可以被**多行 atom_site 数据**描述，每一行都有自己的 `_atom_site_occupancy` 值。一行占位 < 1 表示"这个原子只有占位率这么大的概率在这里"。分两种情况：

**情况 A —— 混合占位**（同一坐标上各行占位之和 = 1.0）：描述*占位无序* —— 位置始终被占据，只是不同物种按概率出现。衍射实验测到的是空间/时间平均结构。

```
Au1 Au 0.0 0.0 0.0 0.5      ← 50% Au
Cu1 Cu 0.0 0.0 0.0 0.5      ← 50% Cu  (坐标相同，占位和=1.0 → 该位满)
```

**情况 B —— 单一物种部分占位**（只一行，占位 < 1.0）：描述*真正的空位*。剩下的 (1 − 占位) 是真的空着。

```
Li1 Li 0.0 0.0 0.0 0.5      ← 50% Li，50% 是空位（隐含）
```

xtalkit 的 `enumerate` 两种情况都处理：
- **情况 A**（多物种混合，占位和为 1.0）：enumlib 直接枚举各物种放在哪里 —— 不需要补充。（若多物种位点占位和 < 1，xtalkit 会先把差额补成空位物种。）
- **情况 B**（单一物种 + 真空位）：xtalkit 自动追加显式 `DummySpecies("X")`，占位 `1 − 占位`，把它内部变成情况 A（`Li0.5 → Li0.5 + X0.5`），然后 enumlib 枚举 Li 与空位的有序排列。

#### 教学案例：往有序 CIF 里引入无序

从数据库（如 Materials Project）下载的 CIF 通常完全有序——所有位点占位都是 1。直接 `xtalkit enumerate` 会得到 0 个结构，因为没有可枚举的无序。要研究无序，得先改 CIF。两种编辑模式覆盖所有情况：

**编辑 1 —— 单物种部分占位（Li/空位，情况 B）。** 把该位点占位从 `1` 改成一个分数；xtalkit 自动补空位。

```
Li  Li0  8  0.229  0.273  0.295  1        →        Li  Li0  8  0.229  0.273  0.295  0.5
```

**编辑 2 —— 多物种混合占位（Ge/P，情况 A）。** 把单行替换成两行，**坐标完全相同**，占位之和=1（pymatgen 会把同坐标两行合并成一个混合位点）。

```
Ge  Ge4  2  0.5  0.5  0.301  1        →        Ge  Ge4   2  0.5  0.5  0.301  0.5
                                               P   P4a   2  0.5  0.5  0.301  0.5
```

**占位必须可整数化。** 多重度 M 的位点占位 p，`M × p` 必须是整数——这就是 enumlib 实际填充的位点数。M=8 的位点上 p 可取 1/8、1/4、3/8、1/2、…；M=2 的位点只能 1/2（或 1）。`0.5` 只是"在 M=8 和 M=2 上都能整数化"的最简单公共值，**不是唯一选择**。若你想要的 p 在 1× 里无法整数化，就增大 `--max-cell-size`（p=1/3 需 3× 超胞）或换个附近的分数。

**CIF 的占位决定了每个输出结构的成分。** enumlib 只打乱"哪些位点放哪种物种"，总数不变，所以每个输出结构的成分都等于你 CIF 定义的成分。两个推论：

- 引入**反位**无序（如 Ge/P）又要保持计量比，就让 Ge 位点和 P 位点**都**改成混合（各 {Ge:0.5, P:0.5}），Ge、P 总数就不变。
- **Li/空位无序一定会改变 Li 总数**（Li 位点改 0.5 就少了 Li）。想在**化学计量成分**下研究 Li 无序，需要劈裂 Li 位点（比有序模型更多的 Li 位点、各自部分占位）——有序 MP CIF 没有这些，得用你要复现的论文里的 Li 位点模型。

**LGPS 具体例子。** 从有序 Li10GeP2S12 CIF 出发（P4_2mc，Z=2；Li0/Li1 为 M=8，Li2/Li3/Ge4/P5/P6 为 M=2）：

- 在 Ge4 和 P6 上做 Ge/P 反位（都 → {Ge:0.5, P:0.5}）：成分守恒在 Li10GeP2S12，1× 可整数化。
- 把 Li0 改成 0.5：该位点变 4 Li + 4 空位（1× 可整数化），但成分掉到 Li8GeP2S12——只有你研究缺锂体系时才这么用。

用 `xtalkit enumerate <改好的>.cif --max-cell-size 1` 跑。改完后务必重新核算成分，确认与你要复现的论文一致，再信任枚举结果。

#### 当 enumlib 返回 0 个结构时

最常见的原因是母相 CIF 的占位在你选的超胞尺寸下**无法整数化**（见上文）。xtalkit 会给出清晰报错。两种修法：

1. **增大 `--max-cell-size`**，让超胞大到能容纳每个物种的整数个数（如占位 1/3 需要 `--max-cell-size 3`）。
2. **准备一份"干净"母相 CIF**，把占位改成接近真实值的有理数。实验 CIF 常报告别扭的分数（如某 Li 位点占位 0.56 = 14/25）。舍入到附近的简单分数（0.56 → 0.5）就得到可枚举的母相；CIF 其他部分（晶胞、其他原子、空间群）原样不动。把占位舍入到附近的有理计量，是枚举无序晶体时标准且通用的近似做法。

**非整数计量比预检。** 跑 `enum.x` 之前，xtalkit 会检查每个物种在 `[--min-cell-size, --max-cell-size]` 范围内是否存在某个 `cell_size` 使其计数（`多重度 × 占位 × cell_size`）为整数。若不存在，直接拒绝并列出每个违规物种及最近的合法分数，例如：

```
[ERR] Cannot enumerate: non-integer stoichiometry (enumlib would likely run away and exhaust memory).
  species  mult  occ      mult*occ   nearest valid (at cell_size 1)
  Li      16    0.6910   11.0560    -> 0.6875
  Li      8     0.6430   5.1440     -> 0.6250
  ...
```

加这道护栏是因为 **enumlib 遇到非整数计量比不会干净失败** —— 它倾向于一直申请内存直到内核杀进程（可能搞崩 WSL/小内存系统，而不是报错）。用建议的舍入值重建 CIF（`xtalkit build --atom-frac ...`）再重跑即可。`--skip-preflight` 可跳过此检查（危险——仅在你确信 enumlib 能处理你的计量比时使用）。

**已知限制：**

- **非整数计量比**：如上，不是简单分数的占位（如 0.56）在小超胞里无法整数化，会得到 0 个结构。舍入到附近分数，或增大 `--max-cell-size`。
- **大规模多位点枚举即便占位干净也可能吃光内存。** 原子数多、无序位点多的原胞（如 LGPS——P42/nmc，50 原子，Li 在 16h+4d+8f 上无序、外加 Ge/P 混合）搜索空间巨大：`enum.x` 单线程，可能要好几 GB 内存和好几分钟。若返回 0 / 崩溃：(a) 设置 `--timeout` 让失控搜索被干净杀掉；(b) **缩小范围**——在母相 CIF 里把不研究的无序位点固定为单一有序（占位 1），一次只枚举一个无序位点；(c) 保持 `--max-cell-size 1`；`--max-structures` 只限制 `enum.x` 之后的输出阶段。
- **平台说明**：`scripts/build_enumlib.sh` 面向 Linux/macOS（系统 `gfortran`）。Windows 用户可在 [WSL](https://learn.microsoft.com/windows/wsl/) 下编译，或沿用旧的 conda + `m2w64-gcc-fortran` 路径（把 `enum.x`/`makestr.x` 放到环境的 `Library/mingw-w64/bin`）。


---

#### 枚举功能配置

xtalkit 核心（`mark`、`skeleton`、`info`、`fetch`、`build`）**不需要** pymatgen，仅 `enumerate` 需要，且它是一个可选 uv extra —— 核心保持轻量（仅 gemmi + rich）。无需 conda，无需 root。

**第 1 步 —— 安装 `enumerate` extra：**

```bash
uv sync --extra enumerate
```

这会装上 `pymatgen`（>=2024.5）。PyPI 上现代版本的 pymatgen 仍内置 `pymatgen.command_line.enumlib_caller.EnumlibAdaptor`，所以无需钉死旧版本。（旧文档里的 `2023.5.31` 钉版是 conda-forge 专用 workaround；uv 走 PyPI，不受影响。）

**第 2 步 —— 编译 enumlib 二进制（一次性，无需 root）：**

```bash
bash scripts/build_enumlib.sh
```

该脚本克隆 [msg-byu/enumlib](https://github.com/msg-byu/enumlib)（含 `symlib` 子模块），用系统 `gfortran` 编译 `enum.x` 与 `makestr.x`，装到 `~/.local/share/xtalkit/bin/`。需要 `gfortran`、`git`、`make`（如 `sudo apt install gfortran make git`）。可用 `XTALKIT_ENUMLIB_BIN` 覆盖安装位置。

xtalkit 会自动发现这些二进制 —— **无需手动配 PATH**。`xtalkit._env.setup_for_enumlib()` 会在 `enumlib_caller` 导入时的 `which("enum.x")` 之前，把安装目录加到 `PATH` 前面。它也会检查 `$XTALKIT_ENUMLIB_BIN`，并为开发用途回退到仓库内 `enumlib_src/enumlib/src/`。

**第 3 步 —— 验证：**

```bash
uv run xtalkit enumerate tests/fixtures/disordered_binary.cif --max-cell-size 2
# 预期：Au0.5/Cu0.5 → 3 个有序 CIF
uv run pytest tests/test_enumerator.py -v
# 预期：enumerator 测试全部通过
```

在 Windows 上，`xtalkit/_env.py` 还会在运行时应用三个 workaround（Linux/macOS 不执行）：

1. 向 `PATHEXT` 追加 `.X` 和 `.PY`，让 `shutil.which("enum.x")` 能找到二进制
2. 调用 `os.add_dll_directory(env/Library/bin)`，让 scipy 的原生扩展能加载
3. Monkey-patch `shutil.which` 使其只返回绝对路径（Windows 默认会返回 `.\makestr.x`，`subprocess.Popen` 无法启动）

---

### `shry` — 严格 SHRY 枚举工作流

`xtalkit shry` 面向 enumlib 可能耗尽内存的大规模部分占据体系，提供分阶段、可审计的 SHRY 工作流。原有 `enumerate` 命令仍保留。

```bash
xtalkit shry prepare input.cif --out ready.cif --parent-spacegroup 137
xtalkit shry count ready.cif --scaling-matrix 1 1 1 --out count.json
xtalkit shry enum ready.cif --expect-count <COUNT> --out shry_enum --write-cif --write-degeneracy
xtalkit shry verify shry_enum --check-count --check-formula --check-dedup --symprec-list 1e-4 1e-3 1e-2
xtalkit shry postprocess shry_enum --shortest-distance --pair Li Li
```

SHRY 通过外部 CLI 调用。若 SHRY 安装在独立环境中，设置 `XTALKIT_SHRY_CMD=/path/to/shry`。完整流程见 [docs/user/shry-enumeration.md](docs/user/shry-enumeration.md)。

---

## TUI（交互模式）

不带参数启动：

```bash
xtalkit
```

```
╔═══════════════════════════════════════════╗
║      xtalkit · Crystal Wyckoff Toolkit   ║
╠═══════════════════════════════════════════╣
║  [1] Mark CIF    — Mark Wyckoff in a     ║
║                     structure            ║
║  [2] Skeleton    — Generate pure Wyckoff ║
║                     skeleton             ║
║  [3] Query SG    — View space group      ║
║                     information          ║
║  [4] Fetch DB    — Verify database       ║
║                     online               ║
║  [5] Enumerate   — Enumerate ordered     ║
║                     configurations       ║
║  [6] Build       — Build CIF from        ║
║                     refinement params    ║
║  [0] Exit                                ║
╚═══════════════════════════════════════════╝
```

### TUI 工作流示例 —— Mark CIF

```
Input your choice: 1

═══════════════════════════════════════════════
  Mark Wyckoff Positions in CIF

  CIF file path: ./Li6PS5Cl.cif
    ✓ found D:\structures\Li6PS5Cl.cif

  Space group number: 216

    Space Group #216: F-43m
    Available Wyckoff positions: 4a  4b  4c  4d  16e  24f  24g  48h

  Wyckoff positions to mark (comma-separated, or 'all'): 4a,24f

  Mode: [1] Overlay  [2] Replace  > 1

  Output format: [1] cif  [2] xyz  [3] all  > 3

  Tolerance in fractional coords (default 0.5): [Enter for default]

  Element override (e.g. '4a:Xe,16e:Kr') or Enter to skip: [Enter]

  Output base path [default: D:\structures\Li6PS5Cl_WYCK]: [Enter]

  ✓ Saved to: D:\structures\Li6PS5Cl_WYCK.cif,
               D:\structures\Li6PS5Cl_WYCK.xyz

  Press Enter to continue...
```

---

## 虚拟原子系统

xtalkit 使用**稀有/惰性元素**，这些元素几乎不会出现在真实 CIF 中：

| 优先级 | 元素 | Z | VESTA 颜色 |
|--------|------|---|------------|
| 1 | Xe（氙） | 54 | 银色 |
| 2 | Kr（氪） | 36 | 浅灰 |
| 3 | Rn（氡） | 86 | 粉色 |
| 4 | Ar（氩） | 18 | 浅蓝 |
| 5 | Ne（氖） | 10 | 浅绿 |
| 6 | He（氦） | 2 | 白色 |

Wyckoff 字母按字母序排序，按优先级依次分配元素；若空间群 Wyckoff 位置超过 6 个则循环复用。

**CIF 标签格式：** `WYCK_<letter>`（如 `WYCK_4a`、`WYCK_16e`）。

**自定义映射：**

```bash
xtalkit mark file.cif --sg 216 --wyckoff 4a,4c --map 4a:He,4c:Ne
```

---

## Overlay 与 Replace 模式

### Overlay（默认）

```
之前：  Li at (0,0,0)   P at (0.25,0.25,0.25)
之后：  Li at (0,0,0)   P at (0.25,0.25,0.25)
        WYCK_4a (Xe) at (0,0,0)
        WYCK_4c (Kr) at (0.25,0.25,0.25)
        ... 其他请求的 Wyckoff 位置的虚拟原子
```

真实原子保留，虚拟原子叠加其上。在 VESTA 中两者都会显示。

### Replace

```
之前：  Li at (0,0,0)   P at (0.25,0.25,0.25)
之后：  WYCK_4a (Xe) at (0,0,0)    ← Li 被替换（匹配到 4a）
        P at (0.25,0.25,0.25)      ← P 保留（未在容差内匹配 4c）
```

匹配到请求 Wyckoff 位置（在容差范围内）的原子**被替换**为虚拟原子，其他原子保持不变。

---

## 匹配容差

容差（`--tol`，默认 0.5）控制一个原子的坐标与 Wyckoff 位置理论坐标之间允许多大偏差，仍视为占据该位置。

- **较大值**（0.5–1.0）：宽松匹配，适合有轻微坐标偏差的实验结构
- **较小值**（0.01–0.1）：严格匹配，只接受非常接近理想 Wyckoff 位置的原子

容差在**分数坐标空间**中应用（不是埃）。对近似等轴的立方晶胞，分数空间 0.5 ≈ 实空间各轴方向 0.5 × a 的距离。

---

## 工作流配方

### 配方 1：研究 Li₆PS₅Cl（SG 216，F-43m）的 Wyckoff 占位

```bash
# 第 1 步：从 Materials Project 下载 CIF
# （假设你已有该文件）

# 第 2 步：查看 F-43m 有哪些 Wyckoff 位置
xtalkit info --sg 216

# 第 3 步：标记全部 Wyckoff 位置，overlay 模式，所有格式
xtalkit mark Li6PS5Cl.cif --sg 216 --wyckoff all --format cif,xyz

# 第 4 步：用 VESTA 打开 Li6PS5Cl_WYCK.cif
# → 全部 8 个 Wyckoff 位置以彩色虚拟原子显示
# → 真实原子（Li、P、S、Cl）仍然显示
# → 可以在 VESTA 中切换原子显示以做对比
```

### 配方 2：创建 Wyckoff 参考骨架

```bash
# 用真实晶胞参数为 F-43m 生成骨架
xtalkit skeleton --sg 216 --wyckoff all \
    --cell "9.85 9.85 9.85 90 90 90" \
    --format cif

# 用 VESTA 打开 SG216_skeleton.cif
# → 直观看到每个 Wyckoff 位置在晶胞中的位置
# → 没有真实原子 —— 纯参考模板
```

### 配方 3：检查哪些原子占据特定 Wyckoff 位置

```bash
# 只标记你关心的 Wyckoff 位置
xtalkit mark structure.cif --sg 225 --wyckoff 4a,8c --mode replace

# Replace 模式下，4a 和 8c 位置的原子会被替换为虚拟原子
# → 在 VESTA 中立即看到："这些位置上是否有原子？"
```

### 配方 4：批处理多个结构

```bash
# 同一空间群下目录中的所有 .cif 文件
for f in *.cif; do
    xtalkit mark "$f" --sg 216 --wyckoff all -o "${f%.cif}_WYCK"
done
```

---

## 已支持的空间群

已提供**全部 230 个空间群**的 Wyckoff 位置数据。内置数据集（`xtalkit/data/wyckoff.json`）源自国际晶体学表 Vol. A，经 [pyxtal](https://pyxtal.readthedocs.io)（MIT）抽取并用 gemmi 对称操作校验 —— 重新生成见 `scripts/build_wyckoff_db.py`。

`xtalkit fetch` 可确认数据集完整（230/230）。

---

## 开发

```bash
uv sync                     # 安装核心依赖
uv sync --extra enumerate   # 同时启用 `enumerate`（拉入 pymatgen）
uv run pytest               # 运行全部测试（未装 enumerate extra 时跳过 5 个）
```

### 项目结构

```
xtalkit/
├── xtalkit/
│   ├── __init__.py      # 包、版本
│   ├── cli.py           # argparse CLI + 5 个子命令
│   ├── tui.py           # 基于 rich 的交互式 TUI
│   ├── spacegroup.py    # Gemmi 空间群查询
│   ├── matcher.py       # 原子 → Wyckoff 位置匹配
│   ├── marker.py        # 核心：在 CIF 中标记 Wyckoff
│   ├── skeleton.py      # 纯 Wyckoff 骨架生成
│   ├── exporter.py      # .cif / .xyz 写入器
│   ├── enumerator.py    # enumlib 封装（延迟导入 pymatgen）
│   ├── _env.py          # enumlib 二进制发现 + Windows 环境修复
│   └── utils.py         # 共用工具
├── tests/
│   ├── fixtures/
│   │   ├── simple.cif         # 测试 CIF（F-43m，Li + P）
│   │   └── disordered_binary.cif  # Au0.5/Cu0.5，用于 enumerate 测试
│   ├── test_spacegroup.py
│   ├── test_matcher.py
│   ├── test_exporter.py
│   ├── test_marker.py
│   ├── test_skeleton.py
│   ├── test_enumerator.py     # enumlib 集成测试（无 pymatgen 时跳过）
│   ├── test_cli.py
│   ├── test_tui.py
│   └── test_integration.py
├── docs/superpowers/
│   ├── specs/2026-06-20-xtalkit-design.md
│   └── plans/2026-06-20-xtalkit.md
├── scripts/
│   └── build_enumlib.sh    # 从源码编译 enumlib（enum.x、makestr.x）
├── pyproject.toml
└── README.md
```

---

## 依赖

| 包 | 用途 | 是否必需 |
|----|------|----------|
| [gemmi](https://gemmi.readthedocs.io/) | 空间群数据、CIF I/O | 是 |
| [rich](https://rich.readthedocs.io/) | TUI 格式化（表格、面板、颜色） | 是 |
| [pymatgen](https://pymatgen.org/) >=2024.5 | `enumerate` 的 enumlib 封装 | 仅 `enumerate` extra |
| [enumlib](https://github.com/msg-byu/enumlib) | 对称不等价构型枚举（Fortran） | 仅 `enumerate` 需要 |
| pytest（开发） | 测试框架 | 是 |

`enumerate` 子命令延迟导入 pymatgen，因此核心工具箱在没有 pymatgen 时也能正常工作。`uv sync --extra enumerate` + `build_enumlib.sh` 路径见 [枚举功能配置](#枚举功能配置)。
