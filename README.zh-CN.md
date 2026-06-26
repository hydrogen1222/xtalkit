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

共五个子命令：

| 命令 | 功能 |
|------|------|
| `mark` | 在 CIF 文件中标记 Wyckoff 位置 |
| `skeleton` | 生成纯 Wyckoff 骨架（不含真实原子） |
| `info` | 查询某个空间群的 Wyckoff 位置信息 |
| `fetch` | 校验空间群数据库完整性 |
| `enumerate` | 枚举对称不等价的有序构型（需 pymatgen + enumlib） |

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
| `--format cif` | cif | 输出格式：`cif`、`vesta`、`xyz`，或逗号分隔如 `cif,vesta,xyz` |
| `-o base` | `{name}_WYCK` | 输出基础路径（每种格式自动加扩展名） |

**示例：**

```bash
# 在立方 F-43m 结构中标记 4a 和 24f（overlay 模式）
xtalkit mark Li6PS5Cl.cif --sg 216 --wyckoff 4a,24f

# 标记全部 8 个 Wyckoff 位置，输出三种格式
xtalkit mark structure.cif --sg 216 --wyckoff all --format cif,vesta,xyz

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
| VESTA | `{name}_WYCK.vesta` | VESTA 原生 XML，含晶胞 + 原子位点。 |
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
| `--format cif` | cif | `cif`、`vesta`、`xyz`，或 `cif,vesta,xyz` |
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
    --cell "5.5 6.3 8.1 90 108.5 90" --format cif,vesta
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
  4a       4      -4         0,0,0
  4b       4      -4         1/2,1/2,1/2
  4c       4      -4         1/4,1/4,1/4
  4d       4      -4         3/4,3/4,3/4
  16e      16     .3m        x,x,x
  24f      24     2..        1/4,0,0
  24g      24     .2.        1/4,1/4,1/4
  48h      48     1          1/4,1/4,1/4
```

---

### `fetch` — 校验数据库

```
xtalkit fetch
```

校验全部空间群数据是否完整。输出：

```
✓ Space group data intact (230/230 OK)
```

---

### `enumerate` — 枚举有序构型

使用 [pymatgen](https://pymatgen.org/) + [enumlib](https://github.com/msg-byu/enumlib)（Hart-Forcade 算法）枚举一个无序 CIF 的所有对称不等价有序构型。该子命令为**可选功能** —— 需要 `enumerate` uv extra 与源码编译的 enumlib 二进制（见下文 [枚举功能配置](#枚举功能配置)）。

```
xtalkit enumerate <input.cif> [options]
```

**必填参数：**

| 参数 | 说明 |
|------|------|
| `input.cif` | 含部分占位/无序占位的 CIF 文件路径（如 Au0.5/Cu0.5） |

**可选参数：**

| 选项 | 默认值 | 说明 |
|------|--------|------|
| `--min-cell-size N` | 1 | 枚举的最小超胞尺寸 |
| `--max-cell-size N` | 2 | 枚举的最大超胞尺寸 |
| `--symm-prec TOL` | 0.1 | 空间群分析的对称容差 |
| `--vacancy-symbol S` | `X` | 空位的 DummySpecies 符号 |
| `--output-dir DIR` | `{name}_enum/` | 枚举 CIF 的输出目录 |
| `--max-structures N` | 不限 | 输出文件数量上限 |
| `--timeout MIN` | 无 | enumlib 子进程超时（分钟） |
| `--format F` | cif | 输出格式（`cif`；`xyz` 为预留） |

**示例：**

```bash
# 在 1×–2× 超胞中枚举 50/50 Au/Cu 二元体系
xtalkit enumerate AuCu_disordered.cif --max-cell-size 2
# 输出：AuCu_disordered_enum/AuCu_disordered_000.cif, _001.cif, _002.cif

# 只保留前 5 个结构
xtalkit enumerate disordered.cif --max-cell-size 3 --max-structures 5

# 自定义输出目录
xtalkit enumerate parent.cif --output-dir ./runs/exp01/
```

**工作流程：**

1. 通过 `pymatgen.core.Structure.from_file` 读取 CIF（原生支持部分占位 + mmCIF）
2. **转换为原胞**（`Structure.get_primitive_structure()`）。这一步至关重要：enumlib 不会自动找原胞，而 F 心/I 心的常规晶胞候选位点数会翻 2×–4× 倍，足以让 enumlib 的 tree_class 数组溢出（如 F-43m 48h 重复度 48，C(48,24) 溢出；原胞里只有 12 个候选位点，C(12,6) = 924，正是能复现文献结果的数量级）。
3. 对任何部分占位位点（单一或多物种），自动追加显式 `DummySpecies(vacancy_symbol)` 表示空位（如 `Li0.5` → `Li0.5 + X0.5`；`Au0.3/Cu0.3` → `Au0.3 + Cu0.3 + X0.4`）
4. 调用 `pymatgen.command_line.enumlib_caller.EnumlibAdaptor`，它会 shell 调用 Fortran 二进制 `enum.x` 与 `makestr.x`
5. 把每个对称不等价的有序构型写为 `<basename>_<NNN>.cif`（以原胞形式输出）

#### CIF 的分数占位是怎么工作的

在 CIF 中，同一个晶体学坐标可以被**多行 atom_site 数据**描述，每一行都有自己的 `_atom_site_occupancy` 值。一行占位 < 1 表示"这个原子只有占位率这么大的概率在这里"。分两种情况：

**情况 A —— 混合占位**（同一坐标上各行占位之和 = 1.0）：描述*占位无序* —— 位置始终被占据，只是不同物种按概率出现。衍射实验测到的是空间/时间平均结构。

```
Au1 Au 0.0 0.0 0.0 0.5      ← 50% Au
Cu1 Cu 0.0 0.0 0.0 0.5      ← 50% Cu  (坐标相同，占位和=1.0 → 该位满)
```

**情况 B —— 单一物种部分占位**（只一行，占位 < 1.0）：描述*真正的空位*。剩下的 (1 − 占位) 是真的空着。

```
Li1 Li 0.3148 0.018 0.6852 0.56    ← 56% Li，44% 是空位（隐含）
```

xtalkit 的 `enumerate` 两种情况都处理：
- **情况 A**（多物种混合，占位和为 1.0）：enumlib 直接枚举各物种放在哪里 —— 不需要补充。（若多物种位点占位和 < 1，xtalkit 会先把差额补成空位物种。）
- **情况 B**（单一物种 + 真空位）：xtalkit 自动追加显式 `DummySpecies("X")`，占位 `1 − 占位`，把它内部变成情况 A（`Li0.56 → Li0.56 + X0.44`），然后 enumlib 枚举 Li 与空位的有序排列。

#### 为什么要造 `Li6PS5Cl_clean.cif`（vs `EntryWithCollCode418490.cif`）

文献给出的 argyrodite CIF（`EntryWithCollCode418490.cif`）报告 Li 在 48h Wyckoff 位置上的占位是 **0.56**（对应 Li6.72PS5Cl 计量）：

```
Li1 Li1+ 48 h 0.3148(19) 0.018(4) 0.6852(19) 0.104(14) 0.56(6) 0
```

`0.56 = 14/25` 无法在任何实际可行的超胞中整数化（需要 25× 超胞 = 15625 个原胞）。enumlib 直接返回 0 个结构。

"干净"版（`Li6PS5Cl_clean.cif`）只改了这一处占位值，把 0.56 舍入成 0.5（对应 Li6PS5Cl 计量）：

```
Li1 Li1+ 48 h 0.3148(19) 0.018(4) 0.6852(19) 0.104(14) 0.5 0
```

0.5（= 1/2）在 1× 超胞里就有整数 Li 数（原胞中 12 个 Li 位点 × 0.5 = 6 个 Li），enumlib 能跑通，产出 **48 个对称不等价的有序构型** —— 与文献一致。CIF 中其他部分（晶胞、P/S/Cl 坐标、空间群）原样不动。

这种舍入是 argyrodite 文献里的标准做法：真实材料是 Li6.72PS5Cl（Li 位置无序），但计算研究枚举的是 Li6PS5Cl（附近的有理计量），因为它可行。如果你确实需要 Li6.72 真实计量，要么用大得多的超胞，要么换一种枚举思路。

**已知限制：**

- **非整数计量比**：占位 0.56（= 14/25）无法在任何小超胞中整数化。enumlib 将返回 0 个结构；xtalkit 会以清晰错误提示，建议增大 `--max-cell-size` 或改用占位为有理数（如把 0.56 舍入为 0.5）的"干净"母相 CIF。对 Li6.72PS5Cl（Li 占位 0.56），可准备一份 Li6PS5Cl 母相 CIF（Li 占位改为 0.5）—— 这样能复现 argyrodite 文献中 ~48 个 Li 有序构型。
- **平台说明**：`scripts/build_enumlib.sh` 面向 Linux/macOS（系统 `gfortran`）。Windows 用户可在 [WSL](https://learn.microsoft.com/windows/wsl/) 下编译，或沿用旧的 conda + `m2w64-gcc-fortran` 路径（把 `enum.x`/`makestr.x` 放到环境的 `Library/mingw-w64/bin`）。

---

#### 枚举功能配置

xtalkit 核心（`mark`、`skeleton`、`info`、`fetch`）**不需要** pymatgen，仅 `enumerate` 需要，且它是一个可选 uv extra —— 核心保持轻量（仅 gemmi + rich）。无需 conda，无需 root。

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
# 预期：8 passed
```

在 Windows 上，`xtalkit/_env.py` 还会在运行时应用三个 workaround（Linux/macOS 不执行）：

1. 向 `PATHEXT` 追加 `.X` 和 `.PY`，让 `shutil.which("enum.x")` 能找到二进制
2. 调用 `os.add_dll_directory(env/Library/bin)`，让 scipy 的原生扩展能加载
3. Monkey-patch `shutil.which` 使其只返回绝对路径（Windows 默认会返回 `.\makestr.x`，`subprocess.Popen` 无法启动）

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

  Output format: [1] cif  [2] vesta  [3] xyz  [4] all  > 4

  Tolerance in Å (default 0.5): [Enter for default]

  Element override (e.g. '4a:Xe,16e:Kr') or Enter to skip: [Enter]

  Output base path [default: D:\structures\Li6PS5Cl_WYCK]: [Enter]

  ✓ Saved to: D:\structures\Li6PS5Cl_WYCK.cif,
               D:\structures\Li6PS5Cl_WYCK.vesta,
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

# 第 3 步：标记全部 Wyckoff 位置，overlay 模式，三种格式
xtalkit mark Li6PS5Cl.cif --sg 216 --wyckoff all --format cif,vesta,xyz

# 第 4 步：用 VESTA 打开 Li6PS5Cl_WYCK.vesta
# → 全部 8 个 Wyckoff 位置以彩色虚拟原子显示
# → 真实原子（Li、P、S、Cl）仍然显示
# → 可以在 VESTA 中切换原子显示以做对比
```

### 配方 2：创建 Wyckoff 参考骨架

```bash
# 用真实晶胞参数为 F-43m 生成骨架
xtalkit skeleton --sg 216 --wyckoff all \
    --cell "9.85 9.85 9.85 90 90 90" \
    --format vesta

# 用 VESTA 打开 SG216_skeleton.vesta
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

目前已提供 38 个空间群的 Wyckoff 位置数据：

| 范围 | 晶系 | 数量 |
|------|------|------|
| 1–2 | 三斜 | 2 |
| 195–230 | 立方 | 36 |

未支持的空间群会抛出 `NotImplementedError` 并附清晰提示。扩展到全部 230 个空间群已在计划中。

**最常见的电池材料（立方空间群）已完全支持。**

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
│   ├── exporter.py      # .cif / .vesta / .xyz 写入器
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
