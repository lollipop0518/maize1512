# 计算样本间遗传距离

**HapMap (hmp.txt) → 0/1/2 矩阵 + QC (snpReady) → 全样本遗传距离矩阵 → 最近20 / 最远20 的 Excel → 热图**

------

## 一、使用说明

1. **0/1/2 编码 & 质控：用 snpReady::raw.data**

   - `raw.data()` 能从碱基编码（A/C/G/T 或 A/B）自动转换成 0/1/2 矩阵，并同时做：
     - 个体缺失率过滤（`sweep.sample`）
     - SNP 缺失率过滤（`call.rate`）
     - MAF 过滤（`maf`）
     - 缺失值填补（多种方法：wright / mean / knni）
   - `outfile = "012"` 时：AA→0，Aa→1，aa→2，正好是我们后续做距离计算最方便的形式。

2. **遗传距离的定义：基于等位基因差异的 Manhattan 距离**

   - 对 0/1/2 编码，$|g_{ik}-g_{jk}|$ 正好是在位点 $k$ 上两个材料的**等位基因数差**（0/1/2 个等位基因不同）。

   - 令 $L$ 为 SNP 数，则：
     $$
     d_{ij}=\frac{1}{2L}\sum_{k=1}^L |g_{ik}-g_{jk}|
     $$
     这就是**所有 SNP 上“不同等位基因数”占总等位基因数的比例**，自然落在 $[0,1]$ 区间，直观又易解释。

   - 在 R 里，用 `dist(geno012, method = "manhattan")` 再除以 `2 * L` 就能一次性算完，效率高、实现简单。

3. **Excel 导出：openxlsx**

   - `openxlsx::write.xlsx()` 直接把 data.frame 写成 `.xlsx`，不依赖 Java，速度和兼容性都很好。

4. **热图：ComplexHeatmap + circlize**

   - 对比 base `heatmap()` / `pheatmap()`，`ComplexHeatmap` 在大矩阵布局、图例、字体、导出 PDF 等方面更细致，顶刊常用。
   - 结合 `circlize::colorRamp2()` 自定义连续色标（如蓝–白–红），可把距离梯度表现得更细腻。