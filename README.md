# BSAtool

This project implements two bulked-segregant analysis (BSA) strategies.

1. **Conventional high- vs. low-bulk comparison.**
2. **Recessive-trait mode.** For single- or two-gene recessive traits, the “low bulk” can be replaced by the average parental SNP-index, calculated as (SNP-index\_Parent1 + SNP-index\_Parent2)/2.

To activate the **recessive-trait mode**, simply replace

```python
LOWPOOL_MODE = "bulk"
```

with

```python
LOWPOOL_MODE = "parents"   # or the equivalent keyword used in the script, e.g. "parents"
```


# BSAtool中文说明

本项目中的BSA支持两种分析策略方法一种是标准的高池和低池，对于单基因或双基因隐性性状混池，还可以使用 (SNP-index\_Parent1 + SNP-index\_Parent2)/2作为低池。
将代码中的

```python
LOWPOOL_MODE = "bulk"
```

改为

```python
LOWPOOL_MODE = "parent_avg"   # 或脚本中定义的等价关键字，例如 "parents"
```

即可启用 **Recessive-trait mode**，让脚本把
$(\text{SNP-index}_{\text{Parent 1}} + \text{SNP-index}_{\text{Parent 2}}) / 2$
作为低池进行比较。



The fourth-power Euclidean-distance statistic (ED⁴) was computed with the custom script **ED4.py**.

