#!/usr/bin/env python3
# ed4_qtl_bsa.py  ——  ED⁴-based QTL detection (parents-vs-bulk 双模式)

"""
1. 读 VCF → 计算 ΔSNP → ED → ED⁴
2. Monte-Carlo / 正态近似给 ΔSNP CI，并取 ED⁴ α 分位作阈值
3. Tricube 平滑，再输出 per-chrom / genome 图表

低池两种模式：
    • LOWPOOL_MODE="parents"  → 低池 = PARENT1+PARENT2 均值   (默认，与旧脚本一致)
    • LOWPOOL_MODE="bulk"     → 低池 = LOW_POOL 指定样本
"""

# ────────────────── 基础库 ──────────────────
import pysam, pandas as pd, numpy as np, matplotlib.pyplot as plt, csv, re, sys
from pathlib import Path
from tqdm import tqdm
import matplotlib as mpl
from math import sqrt
mpl.rcParams["pdf.fonttype"] = mpl.rcParams["ps.fonttype"] = 42
mpl.rcParams["axes.unicode_minus"] = False

# ═════ 用户参数 ═════
VCF_FILE   = "filtered_snps.dp20_all.vcf.gz"

HIGH_POOL  = "ba"            # 高表型池
PARENT1    = "s2"
PARENT2    = "s3"
LOW_POOL     = "za"         # 若 LOWPOOL_MODE="bulk" 才会使用
LOWPOOL_MODE = "bulk"      # "parents" 或 "bulk"
WINDOW_BP  = 1_000_000
STEP_BP    = 100_000
POINT_SIZE = 6
RAW_FRAC   = 0.3
PALETTE    = "tab10"

N_PERM      = 1000
MAX_DP_PERM = 500
ED_ALPHA    = 0.999           # 置换 α 分位阈值
# ═══════════════════

OUT_PREFIX = Path(VCF_FILE).with_suffix("")
OUT_FOLDER = Path("result/ED4")
OUT_FOLDER.mkdir(parents=True, exist_ok=True)

# ── 染色体工具 ──
_ignore = re.compile(r"arahy\.Tifrunner\.gnm2\.scaffold")
def ignore_chr(c): return bool(_ignore.search(c))
def clean_chr(ch):
    ch = ch.replace("arahy.Tifrunner.gnm2.Arahy.","Chr").replace("arahy.YZ9102.","")
    return f"Chr{int(ch[3:]):02d}" if ch.startswith("Chr") and ch[3:].isdigit() else ch

# ── smoothing ──
def tricube_smooth(x, y, win):
    h = win/2
    sm = np.empty_like(y,float)
    for i in range(len(y)):
        w = (1-(np.abs(x-x[i])/h)**3)**3
        w[np.abs(x-x[i])>h]=0
        sm[i] = (w*y).sum()/w.sum() if w.sum() else np.nan
    return sm

# ═════ STEP 1. 遍历 VCF ═════
vcf = pysam.VariantFile(VCF_FILE,"r")
samples_needed = [HIGH_POOL, PARENT1, PARENT2] + ([LOW_POOL] if LOWPOOL_MODE=="bulk" else [])
for s in samples_needed:
    if s not in vcf.header.samples:
        sys.exit(f"[ERROR] sample '{s}' not found in VCF header")

thr_dict = {}
raw_csv = OUT_FOLDER/f"{OUT_PREFIX.name}_ED4_raw.csv"
with open(raw_csv,"w",newline="") as fh:
    wr=csv.writer(fh)
    low_desc = "mean(P1,P2)" if LOWPOOL_MODE=="parents" else f"bulk({LOW_POOL})"
    wr.writerow(["CHROM","POS","high_idx",f"low_idx({low_desc})",
                 "deltaSNP","ED4","CI95","CI99","ED4_thr"])

    for chrom in tqdm([c for c in vcf.header.contigs if not ignore_chr(c)], desc="VCF pass"):
        for rec in vcf.fetch(chrom):
            if len(rec.alts)!=1: continue
            try:
                ad_h  = rec.samples[HIGH_POOL]["AD"]
                ad_p1 = rec.samples[PARENT1]["AD"]
                ad_p2 = rec.samples[PARENT2]["AD"]
                if LOWPOOL_MODE=="bulk":
                    ad_l = rec.samples[LOW_POOL]["AD"]
            except KeyError: continue
            if None in (ad_h,ad_p1,ad_p2) or (LOWPOOL_MODE=="bulk" and ad_l is None): continue

            rh,ah = ad_h; r1,a1 = ad_p1; r2,a2 = ad_p2
            dh,d1,d2 = rh+ah, r1+a1, r2+a2
            if LOWPOOL_MODE=="bulk":
                rl,al = ad_l
                dl = rl+al
                if min(dh,dl)==0: continue
                high_idx = ah/dh
                low_idx  = al/dl
                delta    = high_idx - low_idx
                total_alt = ah+al
                total_dp  = dh+dl
                var_low   = (total_alt/total_dp)*(1-total_alt/total_dp)/dl
            else:  # parents mode
                if min(dh,d1,d2)==0: continue
                high_idx = ah/dh
                low_idx  = (a1/d1 + a2/d2)/2
                delta    = high_idx - low_idx
                total_alt = ah+a1+a2
                total_dp  = dh+d1+d2
                var_low   = (total_alt/total_dp)*(1-total_alt/total_dp)*(1/(4*d1)+1/(4*d2))

            ed  = sqrt(2)*abs(delta)
            ed4 = ed**4
            p = total_alt/total_dp
            var_high = p*(1-p)/dh
            var_total = var_high + var_low
            se = sqrt(var_total)

            if total_dp > MAX_DP_PERM:
                ci95,ci99 = 1.96*se, 2.58*se
                thr = (sqrt(2)*ci95)**4
            else:
                h_perm  = np.random.binomial(dh, p, N_PERM)
                if LOWPOOL_MODE=="bulk":
                    l_perm  = np.random.binomial(dl, p, N_PERM)
                    delta_perm = h_perm/dh - l_perm/dl
                else:
                    p1_perm = np.random.binomial(d1, p, N_PERM)
                    p2_perm = np.random.binomial(d2, p, N_PERM)
                    delta_perm = h_perm/dh - (p1_perm/d1 + p2_perm/d2)/2
                ci95 = np.quantile(np.abs(delta_perm),0.975)
                ci99 = np.quantile(np.abs(delta_perm),0.995)
                thr  = np.quantile((sqrt(2)*np.abs(delta_perm))**4, ED_ALPHA)

            if chrom not in thr_dict:
                thr_dict[chrom] = thr

            wr.writerow([rec.chrom,rec.pos,
                         round(high_idx,5),round(low_idx,5),
                         round(delta,5),round(ed4,5),
                         round(ci95,6),round(ci99,6),
                         round(thr,6)])

vcf.close()

# ═════ STEP 2. 平滑 & 单染色体图 ═════
df=pd.read_csv(raw_csv)
df=df[~df["CHROM"].apply(ignore_chr)]; df["POS"]=df["POS"].astype(int)
smooth_ls=[]; chr_len={}
for ch in sorted(df["CHROM"].unique()):
    sub=df[df["CHROM"]==ch].sort_values("POS").copy()
    chr_len[ch]=sub.POS.max()
    sub["ED4_smoothed"]=tricube_smooth(sub.POS.values, sub.ED4.values, WINDOW_BP)
    sub.to_csv(OUT_FOLDER/f"{OUT_PREFIX.name}_{clean_chr(ch)}_ED4.csv",index=False)
    smooth_ls.append(sub)

    # PDF
    fig,ax=plt.subplots(figsize=(4,6))
    ax.plot(sub.POS,sub.ED4_smoothed,color="tab:red",lw=1.6,label="ED⁴ (smoothed)")
    ax.axhline(thr_dict[ch],color="red",ls="--",label=f"ED⁴ {ED_ALPHA:.1%} thr")
    ax.set(title=f"{clean_chr(ch)}  ED⁴", xlabel="Pos (bp)", ylabel="ED⁴")
    ax.legend(fontsize=8); ax.grid(alpha=0.3); fig.tight_layout()
    plt.savefig(OUT_FOLDER/f"{OUT_PREFIX.name}_{clean_chr(ch)}_ED4.pdf",dpi=300)
    plt.close()

# ═════ STEP 3. 整基因组 ═════
gen=pd.concat(smooth_ls,ignore_index=True)
offset=0;offd={}
for c in sorted(chr_len): offd[c],offset=offset,offset+chr_len[c]+STEP_BP
gen["cumPOS"]=gen.apply(lambda r:r.POS+offd[r.CHROM],axis=1)
cmap=plt.get_cmap(PALETTE); ncol=cmap.N
chr2col={c:cmap(i%ncol) for i,c in enumerate(sorted(chr_len))}
fig,ax=plt.subplots(figsize=(20,6))
for ch in sorted(chr_len):
    sub=gen[gen.CHROM==ch]
    ax.plot(sub.cumPOS,sub.ED4_smoothed,color=chr2col[ch],lw=1.3)
    ax.hlines(thr_dict[ch], offd[ch], offd[ch]+chr_len[ch],
              color="green",ls="--",lw=1)
for ch in sorted(chr_len):
    ax.axvline(offd[ch],color="black",lw=0.3)
ax.set_xticks([offd[c]+chr_len[c]/2 for c in sorted(chr_len)])
ax.set_xticklabels([clean_chr(c) for c in sorted(chr_len)],fontsize=8)
ax.set_xlabel("Chromosomes"); ax.set_ylabel("ED⁴")
ax.set_title(f"Genome-wide ED⁴ (low-pool: {low_desc})")
ax.grid(alpha=0.3); plt.tight_layout()
plt.savefig(OUT_FOLDER/f"{OUT_PREFIX.name}_ED4_allChr.png",dpi=300)
plt.close()

print(f"\n✅ ED⁴ pipeline finished in '{LOWPOOL_MODE}' mode "
      f"(low pool = {low_desc}). Results → {OUT_FOLDER}")
