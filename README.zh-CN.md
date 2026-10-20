# xtalkit

[English README](README.md) · [SHRY 专题手册](docs/user/shry-enumeration.zh-CN.md)

**一句话：把论文里的精修表变成 CIF，再把部分占据的 CIF 变成一批可以拿去算的有序结构。**

```text
精修表 ──build──▶ 无序 CIF ──enumerate / shry──▶ 有序结构 ──ewald──▶ 挑出低能构型 ──▶ DFT
```

顺带还有几个看 Wyckoff 位置的小工具。它不做精修、不做弛豫，也不替你决定氧化态和超胞大小；
遇到这类物理选择时，它会让你自己给，或者直接报错，不会悄悄猜一个。

## 命令一览

| 我想… | 命令 |
|---|---|
| 用精修数据建一个 CIF | `build` |
| 把部分占据结构枚举成有序结构（中小体系） | `enumerate` |
| 同上，但组合数很大，需要先计数再生成 | `shry` |
| 按静电能给一批有序结构排序 | `ewald` |
| 查一个空间群有哪些 Wyckoff 位置 | `info` |
| 在已有 CIF 上标出某些 Wyckoff 位置（只看不算） | `mark` |
| 不需要真实结构，只画出空间群的位置骨架 | `skeleton` |
| 检查内置的 230 个空间群数据是否完整 | `fetch` |

任何命令加 `--help` 都能看到完整参数，例如 `uv run xtalkit build --help`。
直接运行 `uv run xtalkit`（不带子命令）会进入交互式界面，适合记不住参数的时候。

## 安装

```bash
git clone https://github.com/hydrogen1222/xtalkit.git
cd xtalkit
uv sync
uv run xtalkit --version
```

这样就能用 `build`、`info`、`mark`、`skeleton`、`fetch`。
`enumerate`、`shry`、`ewald` 需要额外组件，见后文[可选组件](#可选组件)。

---

## 主线：从精修数据建模

下面以 Li₆PS₅Cl（argyrodite，F-43m，216 号）为例，这类结构里既有半占据的 Li，也有 Cl/S 混占位。

### 第 1 步：把精修表抄成 JSON

推荐写成 JSON 文件，而不是一长串命令行参数：好改、好复查，也能放进课题目录留档。
仓库里有一个可以直接复制的模板：[examples/argyrodite_template.json](examples/argyrodite_template.json)。

```json
{
  "sg": 216,
  "cell": {"a": 9.86, "b": 9.86, "c": 9.86, "alpha": 90, "beta": 90, "gamma": 90},
  "atoms": [
    {"element": "Li", "xyz": [0.182, 0.182, 0.022], "occ": 0.5},
    {"element": "P",  "xyz": [0.5,   0.5,   0.5]},
    {"element": "S",  "xyz": [0.62,  0.62,  0.62]},
    {"element": "Cl", "xyz": [0,     0,     0],    "occ": 0.5},
    {"element": "S",  "xyz": [0,     0,     0],    "occ": 0.5},
    {"element": "S",  "xyz": [0.75,  0.75,  0.75], "occ": 0.5},
    {"element": "Cl", "xyz": [0.75,  0.75,  0.75], "occ": 0.5}
  ]
}
```

> 上面的数值只是示意，用之前请换成你所引用论文的精修值。

精修表和 JSON 字段的对应关系：

| 精修表里的列 | JSON 里写成 | 说明 |
|---|---|---|
| 空间群 | `"sg"` | 只写编号。注意论文用的 setting（原点、轴的选取）要和标准 setting 一致 |
| a, b, c, α, β, γ | `"cell"` | 六个都要写，立方晶系也一样 |
| 原子 / 元素 | `"element"` | 写元素符号，例如 `Li`。不要写 `Li1`、`Li+` |
| x, y, z | `"xyz"` | 直接抄分数坐标即可，不需要换算成标准代表点，程序会自己识别属于哪个 Wyckoff 位置 |
| Occ / SOF | `"occ"` | 0–1 之间，省略时默认为 1 |
| Wyckoff | 不用写 | 程序会自动识别并打印出来，供你和论文核对 |
| Uiso / Biso | 不用写 | 不会写进 CIF |

**混占位**：同一个位置写两行、坐标相同，例如上面的 4a 位就是 `Cl 0.5` 加 `S 0.5`。
同一位置上所有原子的占据之和不能超过 1。

如果论文只给了 Wyckoff 标签和自由参数，也可以这样写：
`{"element": "S", "wyckoff": "16e", "free": [0.62]}`。两种写法可以在同一个文件里混用。

### 第 2 步：运行

```bash
uv run xtalkit build --spec examples/argyrodite_template.json -o LPSC
```

输出大致如下：

```text
[info] SG #216 (F-43m), cubic — detected Wyckoff sites:
        Li  at (0.1820,0.1820,0.0220) occ=0.5   -> 48h (mult 48)
        P   at (0.5000,0.5000,0.5000) occ=1.0   -> 4b (mult 4)
        S   at (0.6200,0.6200,0.6200) occ=1.0   -> 16e (mult 16)
        Cl  at (0.0000,0.0000,0.0000) occ=0.5   -> 4a (mult 4)
        ...
[OK] SG #216 (F-43m), crystal system: cubic, formula: Li6PS5Cl
     Saved to: LPSC.cif
```

### 第 3 步：检查（这一步别跳过）

1. **Wyckoff 是否和论文一致。** 如果识别出的位置和论文不同，或者显示 `?`，
   多半是坐标抄错了，或者论文用的是非标准 setting。
2. **化学式是否正确。** `formula` 是按"多重度 × 占据"算出来的。
   如果你预期 Li₆PS₅Cl，结果却出现小数或比例不对，就回头检查 occ。
   程序**不检查电中性**，这一点要靠你自己。
   化学式里元素的顺序与 JSON 中各行首次出现的顺序一致，想要 `Li6PS5Cl` 就把 Li 写在第一行。
3. 用 VESTA 打开 `LPSC.cif`，看一眼是否合理。

只想临时试一下的话，也可以不写 JSON，直接在命令行里给：

```bash
uv run xtalkit build --sg 216 --cell "9.86 9.86 9.86 90 90 90" \
  --atom-frac "P 0.5 0.5 0.5" \
  --atom-frac "Cl 0 0 0 0.5" --atom-frac "S 0 0 0 0.5" \
  -o test
# 格式：--atom-frac "元素 x y z [占据]"；另有 --atom "元素 Wyckoff [自由参数...] [占据]"，两者可以混用
```

### build 常见报错

| 报错里出现 | 原因和处理 |
|---|---|
| `total occupancy ... exceeds 1.0` | 同一位置占据之和大于 1。精修值四舍五入（如 0.385 + 0.616）就会触发，手动归一化即可 |
| `unit cell is inconsistent` | 晶胞和晶系不匹配，比如立方晶系却 a ≠ c，或者角度不是 90° |
| `Invalid Wyckoff label` | 只在使用 `wyckoff` 写法时出现，标签要写全，例如 `16e`，而不是 `e` |
| `Unknown element` | 元素符号写成了 `Li1`、`li`、`Li+` 等形式 |
| 识别出的位置显示 `?` | 坐标对不上任何 Wyckoff 位置。检查是否抄错，或者空间群 / setting 是否用错 |

---

## 接下来：把无序结构变成有序结构

`build` 得到的是平均结构（带分数占据），DFT 不能直接算。需要选一个超胞，把分数占据"实例化"成整数个原子。

**体系小、组合数不大时，用 `enumerate`**（需要 enumlib）：

```bash
# 先用最小范围试水
uv run xtalkit enumerate LPSC.cif \
  --min-cell-size 1 --max-cell-size 1 \
  --timeout 10 --max-structures 20 --output-dir LPSC_enum
```

- `--min/max-cell-size` 是相对**原胞**的体积倍数，不是沿三个轴分别放大的倍数。
- 如果在当前超胞下占据数凑不成整数，程序会在运行前直接拒绝，并给出最接近的可行占据。
  不要用 `--skip-preflight` 强行绕过：输入有误时 enumlib 可能耗尽内存。

**组合数很大时，用 `shry`**，按 `prepare → count → enum → verify` 分步进行，每一步都先看结果再继续：

```bash
uv run xtalkit shry prepare LPSC.cif --out ready.cif \
  --parent-spacegroup 216 --target-formula Li6PS5Cl --scaling-matrix 1 1 1
uv run xtalkit shry count ready.cif --scaling-matrix 1 1 1 --out count.json
# 看 count.json 里的 count_only_result，确认数量能接受后：
uv run xtalkit shry enum ready.cif --scaling-matrix 1 1 1 --expect-count <数量> \
  --target-formula Li6PS5Cl --write-cif --write-poscar --out SHRY_run
uv run xtalkit shry verify SHRY_run --target-formula Li6PS5Cl
```

完整说明（包括可选的 postprocess）见 [SHRY 专题手册](docs/user/shry-enumeration.zh-CN.md)。

> **argyrodite 特别注意**：48h 的 Li 位成对出现，同一对的两个位置相距不到 1 Å，
> 物理上不可能同时有 Li。枚举工具并不知道这一点，会生成很多这种"两个位置都被占"的构型。
> 枚举之后最好按最短 Li–Li 距离过滤一遍，再进入下一步。

## 最后：用 Ewald 能挑出低能构型

```bash
uv run xtalkit ewald LPSC_enum/ \
  --charges Li:1 P:5 S:-2 Cl:-1 \
  --top-n 20 --group --group-dir ewald_results
```

- 结果按静电能排序，完整排名写进 CSV；`--group` 会把前 N 个**复制**到 `selected/`（只有加 `--move` 才会移动原文件）。
- 只有成分相同、原子数相同的结构才能放在一起排；必须电中性。
- Ewald 能只是初筛手段，不能当作 DFT 能量。

---

## 其他小工具

```bash
uv run xtalkit info --sg 216                                  # 列出 Wyckoff 位置
uv run xtalkit mark input.cif --sg 216 --wyckoff 4a,48h -o m  # 在 CIF 上标出位置
uv run xtalkit skeleton --sg 216 --wyckoff all \
  --cell "9.86 9.86 9.86 90 90 90" -o skel                    # 画位置骨架
uv run xtalkit fetch                                          # 自检内置数据库
```

`mark` 和 `skeleton` 的输出带有 `visualization_only` 标记，**只用于观察**，不能拿去做计算。
`mark` 会核对 `--sg` 和 CIF 里声明的空间群是否一致；不一致时请先查清 setting，不要为了让命令通过而改编号。

## 可选组件

### enumlib（`enumerate` 需要）和 pymatgen（`ewald` 需要）

```bash
sudo apt install gfortran make git      # 编译依赖
uv sync --extra enumerate
bash scripts/build_enumlib.sh           # 安装到 ~/.local/share/xtalkit/bin/，无需改 PATH
```

自测：

```bash
uv run xtalkit enumerate tests/fixtures/disordered_binary.cif \
  --max-cell-size 2 --output-dir enum_smoke     # 应得到 3 个有序 CIF
```

构建脚本默认固定到已审计的 enumlib 提交（`67217859c1cf…`），默认单进程编译。
确有需要时，可以用 `ENUMLIB_REF=<提交>` 换源码版本，用 `XTALKIT_BUILD_JOBS=<N>` 增加编译并发。

### SHRY（`shry` 需要）

SHRY 1.1.8 依赖旧版 pymatgen，**必须单独安装**，不要用 `uv pip install shry` 装进本项目：

```bash
uv python install 3.11
uv tool install --python 3.11 'shry==1.1.8'
shry --help
# 若集群上找不到 shry：export XTALKIT_SHRY_CMD="$(command -v shry)"
```

### tblite（可选）

xtalkit 只负责生成 tblite 作业目录，不自带可执行文件。推荐装在独立的 conda-forge 前缀里
（[安装说明](https://tblite.readthedocs.io/en/stable/installation.html)）：

```bash
MAMBA_ROOT_PREFIX="$HOME/.local/share/micromamba" \
  micromamba create --yes --prefix "$HOME/.local/share/xtalkit/tblite-0.7.0" \
  --channel conda-forge --strict-channel-priority 'tblite=0.7.0'
ln -sfn "$HOME/.local/share/xtalkit/tblite-0.7.0/bin/tblite" "$HOME/.local/bin/tblite"
```

提交作业前要按实际体系确认总电荷和自旋，不要默认 `--charge 0 --spin 0` 总是对的。

## 输出文件

| 命令 | 输出 |
|---|---|
| `build` | 物理 CIF（带对称操作、占据和化学式），可选 XYZ |
| `enumerate` | 默认写到 `<输入名>_enum/`，每个文件是一种有序构型 |
| `shry enum` | `clean_cif/` 是去掉空位后的结构；`manifest.json` 记录本次运行；`raw_shry/` 是原始结果 |
| `ewald` | 排名 CSV；加 `--group` 时另有 `selected/` 和 `rest/` |
| `mark`、`skeleton` | 只用于可视化的 CIF/XYZ |

同名文件会被覆盖，重要任务请显式指定 `-o` 或输出目录。`shry enum` 的输出目录必须是新目录或空目录。

## 常见问题

- **`VIRTUAL_ENV=... does not match`**：当前 shell 激活了别的虚拟环境，先运行 `deactivate`，再用 `uv run`。
- **找不到 `enum.x`**：重新运行 `bash scripts/build_enumlib.sh`。
- **`non-integer stoichiometry`**：在当前超胞大小下，占据数凑不成整数个原子。加大 `--max-cell-size`，或者回头核对占据。
- **空间群对不上**：核对 CIF 里的 H-M 符号、编号、setting 和对称操作是否自洽。程序不会替你猜。
- **SHRY 数量巨大**：回到 `shry count`，检查是否误用了 F/I 心的常规胞、缩放矩阵是否设对。
- **Ewald 报不电中性**：检查氧化态和化学式。程序不会替你补背景电荷。

## 开发

```bash
uv sync --extra enumerate
uv run pytest
uvx ruff check xtalkit tests
```
