# GBLUP模型

## 简介

​        传统意义上的 kinship 矩阵（亲缘系数矩阵）是根据个体之间的血缘关系或遗传概率定义的，而 GRM 则是用基因型数据估计的总体遗传相似性矩阵。在很多群体遗传学、基因组选择或 GWAS（基因组关联研究）中，常用 GRM 作为 kinship 结构的替代或输入，以调整样本间的相关性。

**生成基因组关系矩阵（GRM）**
使用 `G.matrix()` 函数可以根据 SNP 标记数据计算不同类型的基因组关系矩阵（如 VanRaden 方法等）。这些矩阵反映了样本之间的遗传相似性，与亲缘关系矩阵的目的类似（用于遗传分析、基因组选择等）。

```
library(snpReady)
data(maize.hyb)
x <- G.matrix(maize.hyb, method="VanRaden", format="wide")
A <- x$Ga  # additive genomic relationship matrix
D <- x$Gd  # dominance genomic relationship matrix (若method支持)
```

在 **snpReady::G.matrix()** 的输出中：

```
A <- x$Ga
D <- x$Gd
```

这两个矩阵分别代表 **不同遗传效应层面的基因组关系矩阵** 👇

------

## 1️⃣ `Ga`：加性基因组关系矩阵（Additive GRM）

### 含义

**`Ga` 表示加性遗传效应（additive genetic effects）的基因组关系矩阵**。

- 衡量的是个体之间 **等位基因替换效应** 的相似性
- 对应经典数量遗传学中的 **加性亲缘关系矩阵 A**
- 是 **最常用、最核心** 的 kinship / GRM

### 生物学意义

- 假设每个位点的效应是等位基因 **线性相加**
- 不考虑显性或上位性
- 适合：
  - GBLUP
  - 混合线性模型
  - GWAS 中的 kinship 矫正

### 数学直觉（VanRaden 方法）

$$
G_A = \frac{ZZ^\top}{2\sum p(1-p)}
$$

其中：

- $Z$：中心化的 SNP 基因型矩阵
- $p$：等位基因频率

### 实际用途

```
y = Xb + Zu + e
u ~ N(0, Ga * σ²_A)
```

👉 **如果只需要“kinship 矩阵”，通常就是用 `Ga`**

------

## 2️⃣ `Gd`：显性基因组关系矩阵（Dominance GRM）

### 含义

**`Gd` 表示显性遗传效应（dominance effects）的基因组关系矩阵**。

- 衡量的是个体之间 **杂合型（Aa）效应** 的相似性
- 描述等位基因之间的 **非线性互作（显性）**

### 生物学意义

- 反映 **杂交优势（heterosis）**
- 在杂交作物（玉米、水稻、油菜等）中尤为重要
- 不体现在传统 pedigree A 矩阵中

### 模型中的作用

```
y = Xb + Za + Zd + e
a ~ N(0, Ga * σ²_A)
d ~ N(0, Gd * σ²_D)
```

### 什么时候需要？

- 研究 **杂交优势**
- 需要同时估计：
  - 加性方差 σ²_A
  - 显性方差 σ²_D

------

## 3️⃣ Ga vs Gd 总结对比

| 项目             | `Ga`         | `Gd`         |
| ---------------- | ------------ | ------------ |
| 遗传效应类型     | 加性效应     | 显性效应     |
| 是否等同 kinship | ✅ 是（常用） | ❌ 否         |
| 是否常规分析必需 | ✅            | ⭕ 进阶       |
| 是否用于 GBLUP   | ✅            | 仅在 AD 模型 |
| 是否解释杂交优势 | ❌            | ✅            |

------

## 4️⃣ 实务建议

- **GWAS / 群体结构校正 / 常规 GS**

  ```
  K <- x$Ga
  ```

- **杂交群体 / 杂种优势分析**

  ```
  K_A <- x$Ga
  K_D <- x$Gd
  ```

