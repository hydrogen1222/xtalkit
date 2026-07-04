# Plan: `xtalkit enumerate` — streaming + parallel structure generation

**Date:** 2026-07-04
**Status:** proposed

## 目标

降低 `xtalkit enumerate` 的内存占用与磁盘 I/O，并对结构生成阶段（makestr+解析+写文件）做多核并行。`enum.x` 枚举核心保持串行（不按晶胞大小拆分，用户选 B）。

## 瓶颈分析（已确认）

pymatgen `EnumlibAdaptor.run()` 三步：
1. `enum.x`（`_run_multienum`）—— 串行 Fortran 枚举，写 `struct_enum.out`。**单核，不可并行**。
2. `makestr.x struct_enum.out 0 N-1`（`_get_structures`）—— 一次性生成**全部** N 个 `vasp.*` 文件。可按索引区间拆分。
3. `glob("vasp.*")` → 把**全部** vasp 文件解析成 Structure，全部留在 `self.structures` 内存里。可并行 + 内存大头。

当前 `--max-structures` 在第 3 步**之后**才切片 —— 所以它现在不减内存/磁盘。`glob` 返回的是文件名词法序（vasp.0, vasp.1, vasp.10, vasp.11, vasp.2 …），并非数值序，输出顺序其实不稳定。

## 方案（选项 B）

绕过 `adaptor.run()`，自己调 `adaptor._gen_input_file()` + `adaptor._run_multienum()` 拿到 `struct_enum.out` + `num_structs`（这两个方法用 CWD，所以在自己的 scratch 目录里跑）。`_gen_input_file` 会设好 `adaptor.index_species` 与 `adaptor.ordered_sites`，复用之。

### 1. 分批流式 makestr+解析+写 CIF
- `effective_n = min(num_structs, max_structures or num_structs)` —— **生成阶段就限额**（不再先全生成再切）。
- 把 `[0, effective_n)` 切成 `batch_size`（默认 256）的块。每块：
  - 在该块专属临时目录跑 `makestr.x <abs struct_enum.out> <start> <end>`（区间生成，vasp 文件不冲突）。
  - 用 `_parse_vasp_to_structure(data, index_species, ordered_sites)` 解析每个 vasp 文件（**忠实复制** pymatgen `_get_structures` 的单文件逻辑：regex 修正 `scale factor`/负数空格、`Poscar.from_str`、ordered_sites 超胞映射、跳过空位 X、`Structure.from_sites(sorted)`）。
  - 按全局索引写 CIF/XYZ 到 `out_dir`，删掉该块临时目录（释放 vasp 文件 + 内存）。
- 内存上限 ≈ 一个 batch；磁盘只留一个 batch 的 vasp 文件。

### 2. 并行（`--jobs`）
- `jobs=1`（默认）：串行流式（已解决内存/磁盘）。
- `jobs>1`：用 `concurrent.futures.ProcessPoolExecutor` 并行处理各块。worker 函数 `_process_chunk(...)` 模块级（可 pickle）：makestr.x + 解析 + 写文件。Linux fork 继承 PATH（已 `setup_for_enumlib`）+ pymatgen 已导入。
- `jobs=0` = `os.cpu_count()`。
- 注意：每个 worker 进程会占一份 pymatgen（~200MB），jobs 越大内存越大 —— 文档里说明这个权衡。

### 3. 内存盘 scratch（`--scratch-dir`）
- `scratch_dir=None`（默认，系统 temp）或指向 `/dev/shm`（tmpfs，struct_enum.out + vasp 中间件全在内存，避免磁盘 I/O）。
- 注意 `/dev/shm` 通常限 RAM 的 50%，struct_enum.out 可能上 GB —— 文档说明放不下就回退磁盘。

### 4. 确定性
- 输出按**数值索引**命名（`basename_000, basename_001, …`），各块写不相交的索引区间，最终按索引排序。比当前 `glob` 词法序更稳定、可复现。

## 改动文件

- `xtalkit/enumerator.py`：重构 `enumerate_structures` —— 绕过 `run()`，分批流式 + 可选并行；新增模块级 `_parse_vasp_to_structure()` 与 `_process_chunk()` worker；新参数 `jobs`/`batch_size`/`scratch_dir`。保留 `--max-structures`/`--timeout`/`--format` 等。
- `xtalkit/cli.py`：`enumerate` 子命令加 `--jobs`（默认 1）、`--batch-size`（默认 256）、`--scratch-dir`（默认 None）。
- `xtalkit/tui.py`：`_enumerate_workflow` 加 jobs/batch-size/scratch-dir 提示（带默认，可回车跳过）。
- `tests/test_enumerator.py`：① 分批输出与非分批等价（结构集合相同，对比小例的组成/数量）；② `--max-structures` 真正限制生成量（生成数 ≤ 限额）；③ `--jobs=2` 与 `--jobs=1` 结果一致（同结构集合，按索引排序后比较）；④ `_parse_vasp_to_structure` 对照 pymatgen `_get_structures` 同一 vasp 文件解析结果一致。
- `README.md` + `README.zh-CN.md`：文档化 `--jobs`/`--batch-size`/`--scratch-dir`；说明 **enum.x 是串行的**（并行只作用于结构生成阶段；若 enum.x 占大头，并行提速有限）；推荐 `/dev/shm`。

## 风险与缓解

- **复制 pymatgen vasp 解析**：逐行复制 `_get_structures` 单文件逻辑；测试④直接对照 pymatgen 解析同一文件，确保一致。
- **makestr.x 区间语义**：Fortran 版 `makestr.x struct_enum.out start end`（含 end），`.py` 版 `-input struct_enum.out 1 N` —— 按 pymatgen 现有判断（`".py" in MAKESTR_CMD`）分支处理。
- **multiprocessing pickle**：`_process_chunk` 模块级；参数均为可 pickle 类型（字符串/整数/`index_species` 字符串列表/`ordered_sites` 的 `PeriodicSite` 列表可 pickle）。
- **worker 内存**：每 worker 一份 pymatgen；文档提示 jobs 与内存的权衡，默认 1。
- **向后兼容**：默认 `jobs=1`/`batch_size=256`/`scratch_dir=None`，输出结构集合与当前一致（且索引顺序更稳定）。现有测试应通过。

## 不做（本计划外）

- `enum.x` 本身并行（按晶胞大小拆分，选项 C）—— 用户未选；留作后续。
- 改 enumlib Fortran 源码 —— 不动。

## 验证

- 现有 `tests/test_enumerator.py` 全过（行为兼容）。
- 新增等价性/限额/并行一致性/解析对照测试。
- 用用户的 LGPS 无序 CIF 实跑：对比改前改后的结构数量、内存峰值（`/usr/bin/time -v`）、墙钟时间；确认 `--jobs` 与 `/dev/shm` 的效果。
