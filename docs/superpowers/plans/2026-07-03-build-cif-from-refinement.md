# Plan: `xtalkit build` — 从 XRD 精修参数生成 CIF

**Date:** 2026-07-03
**Status:** proposed

## 目标

新增 `build` 子命令：由精修结果（空间群 + 晶胞参数 + 每个原子的 Wyckoff 位点 + 自由坐标值 + 占据率）拼装出标准 CIF。覆盖全部 230 个空间群，运行时仅依赖 gemmi（无新增运行时依赖）。

## 已确认的关键决策

- **输入模型**（用户已确认）：空间群号 + 晶胞 (a,b,c,α,β,γ) + 每原子 (元素, Wyckoff 标签, 自由坐标值, 占据率)。
  - 晶系由空间群自动推导（用于校验晶胞约束，不单独问用户）。
  - 占据率默认 1.0；支持部分占据/混合占据（无序 → 可直接喂给 `enumerate`）。
- **Wyckoff 数据**（用户已确认"内置完整 230 数据集"）：构建期从 pyxtal（MIT）抽取全部 230 个空间群的 Wyckoff 位点，生成内置 JSON。运行时仅 gemmi。**附带收益**：`mark`/`skeleton`/`info`/`fetch` 从 38 个空间群升级到 230 个。
- **接口**（用户已确认）：CLI 参数 + JSON spec + TUI 向导，三者全做。
- **命令名**：`build`（动词，与 mark/skeleton/enumerate 一致）。
- **运行时依赖**：仅 gemmi（核心依赖，已存在）。`build` 是**核心命令**，无需 `--extra`（区别于需要 pymatgen+enumlib 的 `enumerate`）。pyxtal 仅在构建脚本中使用。

## 探查结论（已验证）

- pymatgen `SpaceGroup`、spglib 2.7.0、gemmi 0.7.5 **均不含** Wyckoff 位点坐标表；只能检测已有原子的 Wyckoff 字母。
- pyxtal 1.1.4 的 `Wyckoff_position` 提供：`letter`、`multiplicity`、`get_label()`、`get_dof()`、`get_position_from_free_xyzs(values)`。`get_dof`/`get_position_from_free_xyzs` 在三斜/单斜/六方/立方各晶系均正确（含 `x,2x,z`、`1/3,2/3,z` 等复合表达式）。
- POC 已证明：通过线性回归探测 `get_position_from_free_xyzs` 可还原规范坐标模板串（如 `x,-x,z`、`1/3,2/3,z`、`x,2x,z`），存入 JSON 后运行时用自研简单求值器即可，无需运行时 pyxtal。

## Phase 1 — Wyckoff 数据底座（让全工具支持 230）

1. **`scripts/build_wyckoff_db.py`**（仅构建期；`uv run --with "pyxtal>=1.1" python scripts/build_wyckoff_db.py`）：
   - 对 SG 1–230 的每个 Wyckoff 位点抽取：letter、label、multiplicity、site_symmetry、坐标模板。
   - 模板推导：用 `get_dof()` 得自由参数数 D，对 D 个 slot 各赋不同哨兵值调用 `get_position_from_free_xyzs`，线性回归得每个输出坐标 = Σ cᵢ·slotᵢ + const，规范化为 `x,-x,z`、`1/3,2/3,z`、`x,2x,z`（负分数归一化到 [0,1)；自由参数按物理轴命名）。
   - **逐位点校验**：代入样例值 → gemmi 对称操作展开 → 去重 → 断言 orbit 大小 == multiplicity。
   - **38 个手编 SG 交叉校验**：与现有 `_WYCKOFF_DB` 比对，捕获 setting 差异。
   - `is_standard_setting()` 检查，非标准 setting 标记/跳过。
2. **`xtalkit/data/wyckoff.json`**（提交入库）：`{"<sg>":[{"letter","multiplicity","site_symmetry","coordinates"}, ...]}`。附 pyxtal/MIT + Bilbao/ITA 出处注释。
3. **重构 `xtalkit/spacegroup.py`**：`wyckoff_positions()` 改为加载 `wyckoff.json`（替换 `_WYCKOFF_DB`）。公共 API 不变（`WyckoffInfo`、`wyckoff_positions`、`default_cell_params`、`crystal_system`、`sg_name`）。`NotImplementedError` 仅对真正缺失的 SG 触发（预期无）。
4. 更新 `fetch` 输出（230/230）与 `tests/test_spacegroup.py`（断言 230；跨晶系抽检若干 SG）。
   - **本阶段可独立交付**：mark/skeleton/info 即支持 230。

## Phase 2 — `build` 核心（CLI + 写入器）

5. **`xtalkit/builder.py`**：
   - `free_params(coord_str) -> list[str]`：按出现顺序提取 distinct x/y/z。
   - `eval_coord(coord_str, values) -> (x,y,z)`：线性求值器（处理 `x`、`-x`、`2x`、`1/2-x`、分数、字面量），结果 wrap 到 [0,1)。
   - `build_structure(sg, cell, atoms) -> gemmi.Structure`：每原子校验 letter+元素，`eval_coord` 代入自由值，设 `atom.occ`，把代表原子放入单一 model/chain/residue（仿 `skeleton.py`）。
   - `validate(...)`：晶胞 vs 晶系约束（warn 不 error）、自由参数数 == dof、占据率 ∈ (0,1]、同位点总占据 ≤ 1.0、重复代表原子检测（两用户原子在 SG 下对称等价）。
   - `stoichiometry(structure, sg) -> dict[elem,count]`：复用 `matcher._apply_symmetry` 展开代表原子 + 去重，× 占据率，格式化化学式（如 `Li6PS5Cl`）。
6. **`xtalkit/exporter.py`** 增 `write_structure_cif(structure, path)`：标准 CIF（`_cell_*`、`_symmetry_space_group_name_H-M`、`_symmetry_equiv_pos_as_xyz` loop、`_atom_site_{label,type_symbol,fract_x,fract_y,fract_z,occupancy}`），仅写不对称单元，由读入方展开。增 `write_structure_xyz(structure, sg, path)`：对称展开全晶胞（辅助）。
7. **`xtalkit/cli.py`** 增 `build` 子命令：`--sg`（必）、`--cell "a b c α β γ"`（必）、`--atom` 可重复（`"Li 16e 0.25"` 或 `"Li 16e 0.25 1.0"`：尾部数字按顺序为自由参数，若数量 == dof+1 则末位为占据率）、`--spec path.json`、`--format cif[,xyz]`、`-o/--output`。`cmd_build` 编排 + 打印化学式 + 读回校验。
8. **JSON spec schema**（`--spec`）：`{sg, cell:{a,b,c,alpha,beta,gamma}, atoms:[{element, wyckoff, free:[...], occ}]}`。

## Phase 3 — TUI 向导

9. **`xtalkit/tui.py`** 增 `_build_workflow`：SG → 展示 Wyckoff 菜单（label/mult/site-sym/coords，现 230）→ 晶胞（预填晶系默认值，显示约束）→ 添加原子循环（选元素、选字母、按坐标模板逐个提示自由参数、占据率默认 1.0）→ 实时显示组分表 → 写出。主菜单加 `[6] Build`。
10. 更新 `tests/test_tui.py`。

## Phase 4 — 测试、文档、示例

11. 测试：`tests/test_builder.py`（`free_params`/`eval_coord` 对全 230 模板、构建 NaCl/argyrodite/钙钛矿、化学式、校验报错）、`tests/test_cli.py`（`build` flags + spec）、`tests/test_exporter.py`（结构 CIF 经 gemmi 读回 + VESTA 式展开）、集成测试（build → 读回 → 校验 cell/SG/原子数/化学式）。
12. 文档：README + README.zh-CN 增"Build CIF from refinement"章节 + 示例（如 NaCl：`--sg 225 --cell 5.64 5.64 5.64 90 90 90 --atom "Na 4a" --atom "Cl 4b"`）。更新 `docs/superpowers/specs` 设计文档。

## 风险与缓解

- **模板推导边界情况**（六方 `2x`、`-x`、settings）：构建脚本对**每个**位点用 orbit 大小 == multiplicity 校验；38 个 SG 与手编库交叉比对。问题 SG 优雅降级为 unsupported 并记录。
- **Setting 对齐**（原点选择、单斜轴）：统一用 gemmi 参考设置；构建脚本 `is_standard_setting()` 检查；v1 文档说明仅支持标准 setting。
- **pyxtal 版本漂移**：构建命令 pin `pyxtal>=1.1`；JSON 为权威产物，仅更新时重建。
- **CIF 读回正确性**：集成测试用 gemmi + pymatgen 读回，校验全晶胞原子数与组分。

## 暂不做（v1 范围外）

- 非标准 setting / 原点选择。
- 直接读取精修输出文件（GSAS-II/FullProf/TOPAS）——Q1 中已延后。
- ADP/B-iso、各向异性位移参数输出。
- 磁空间群。
- XYZ 展开输出为次要；若需时间盒，先只做 CIF。
