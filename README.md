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




The fourth-power Euclidean-distance statistic (ED⁴) was computed with the custom script **ED4.py**.

