#!/usr/bin/env python3
# HAT-FLEX Python Code

# --- imports ---
import argparse, gzip, io, re, sys, tarfile, math, json, statistics, time, os, shutil
from collections import Counter, deque, defaultdict
from pathlib import Path
from typing import List, Tuple, Optional, Dict, Iterable, Set

# --- optional pysam (bgzip/tabix) ---
try:
    import pysam  # type: ignore
    HAVE_PYSAM = True
except Exception:
    HAVE_PYSAM = False

# -----------------------------
# Utilities & logging
# -----------------------------
def log(msg: str, v: bool):
    if v:
        print(msg, file=sys.stderr)

def info(msg: str):
    print(msg, file=sys.stderr)

# -----------------------------
# Sex / PAR helpers
# -----------------------------
class ParIntervals:
    """PAR intervals (1-based inclusive); tolerant X/chrX/Y/chrY. Custom BED (0-based half-open) supported."""
    def __init__(self, par_bed: Optional[str]=None, verbose: bool=False):
        self.par = {"X": [], "Y": [], "chrX": [], "chrY": []}
        if par_bed and Path(par_bed).is_file():
            info(f"[PAR] reading pseudoautosomal regions from: {par_bed}")
            self._load_bed(par_bed, verbose)
        else:
            info("[PAR] using built-in hg38 PAR defaults")
            for chrom in ("X","chrX"):
                self.par[chrom].extend([(60001, 2699520), (154931044, 156030895)])
            for chrom in ("Y","chrY"):
                self.par[chrom].extend([(10001, 2649520), (59034050, 59363566)])
        if verbose:
            n = sum(len(v) for v in self.par.values())
            log(f"[PAR] loaded intervals (total entries: {n})", verbose)

    def _load_bed(self, bed_path: str, verbose: bool):
        opener = gzip.open if str(bed_path).endswith(".gz") else open
        with opener(bed_path, "rt") as fh:
            for line in fh:
                if not line.strip() or line.startswith("#"): continue
                p = line.split()
                if len(p) < 3: continue
                chrom, s, e = p[0], p[1], p[2]
                try: s0, e0 = int(s), int(e)
                except ValueError: continue
                if s0 > e0: s0, e0 = e0, s0
                s1, e1 = s0+1, e0
                keys = {chrom}
                keys.add(chrom[3:] if chrom.startswith("chr") else "chr"+chrom)
                for k in keys:
                    self.par.setdefault(k, []).append((s1, e1))

    def in_par(self, chrom: str, pos: int) -> bool:
        L = self.par.get(chrom) or self.par.get(chrom[3:] if chrom.startswith("chr") else "chr"+chrom)
        if not L: return False
        for s,e in L:
            if s <= pos <= e: return True
        return False

def chrom_base(chrom: str) -> str:
    return chrom.replace("chr","")

def child_is_haploid(child_sex: Optional[str], chrom: str, pos: int, par: ParIntervals) -> bool:
    c = chrom_base(chrom)
    if c not in ("X","Y"): return False
    if (child_sex or "").lower().startswith("m"):
        return not par.in_par(chrom, pos)
    return False

# -----------------------------
# VCF helpers
# -----------------------------
def parse_gt(sample_field: str) -> Optional[str]:
    if not sample_field or sample_field == ".": return None
    gt = sample_field.split(":")[0].strip()
    if gt in (".","./.",".|.","." ): return None
    return gt

def parse_gt_alleles(gt: Optional[str]) -> Optional[List[int]]:
    if gt is None or gt in (".","./.",".|."): return None
    if "/" in gt or "|" in gt:
        try: return [int(a) for a in gt.replace("|","/").split("/") if a!="."]
        except ValueError: return None
    try: return [int(gt)] if gt!="." else None
    except ValueError: return None

def is_hom_ref(gt: Optional[str]) -> bool:
    if gt is None: return False
    if gt in ("0","0/0","0|0"): return True
    parts = gt.replace("|","/").split("/")
    return len(parts)==2 and parts.count("0")==2

def is_non_ref(gt: Optional[str]) -> bool:
    if gt is None: return False
    parts = gt.replace("|","/").split("/") if ("/" in gt or "|" in gt) else [gt]
    return any(a not in ("0",".") for a in parts)

def parse_format_sb(sample_field: str, fmt_keys: List[str]) -> Optional[List[int]]:
    if "SB" not in fmt_keys: return None
    i = fmt_keys.index("SB")
    parts = sample_field.split(":")
    if i >= len(parts): return None
    vals = parts[i].replace(",", " ").split()
    try:
        arr = [int(x) for x in vals]
        return arr if len(arr)==4 else None
    except ValueError:
        return None

def alt_on_both_strands(sb: Optional[List[int]], min_each: int=1) -> bool:
    if not sb or len(sb)!=4: return True
    return (sb[2] >= min_each) and (sb[3] >= min_each)

def is_bad_homopolymer_param(ref: str, alt: str, max_run: int) -> bool:
    if max_run <= 1: return False
    runA = "A"*max_run; runT = "T"*max_run
    return (runA in ref or runT in ref or runA in alt or runT in alt)

# Fisher exact right-tail for Strand Bias
def fisher_right_tail(a,b,c,d):
    import math
    n = a+b+c+d
    row1 = a+b; row2=c+d; col1=a+c
    def lC(nn,kk): return math.lgamma(nn+1)-math.lgamma(kk+1)-math.lgamma(nn-kk+1)
    def hyper(x): return math.exp(lC(row1,x)+lC(row2,col1-x)-lC(n,col1))
    x_max = min(col1, row1); p = 0.0
    for x in range(a, x_max+1): p += hyper(x)
    return min(1.0, max(0.0, p))

# -----------------------------
# Allele normalization & keys
# -----------------------------
SYMBOLIC_ALLELES = {"<NON_REF>", "<*>"}

def is_symbolic(alt: str) -> bool:
    """Return True if ALT is symbolic or a breakend notation."""
    return (alt in SYMBOLIC_ALLELES) or alt.startswith("<") or ("]" in alt) or ("[" in alt)

def normalize_allele(pos: int, ref: str, alt: str) -> Tuple[int,str,str]:
    """Trim common prefix/suffix per VCF rules; keeps >=1 base in ref and alt."""
    # Trim common prefix
    while len(ref)>1 and len(alt)>1 and ref[0]==alt[0]:
        ref = ref[1:]; alt = alt[1:]; pos += 1
    # Trim common suffix
    while len(ref)>1 and len(alt)>1 and ref[-1]==alt[-1]:
        ref = ref[:-1]; alt = alt[:-1]
    return pos, ref, alt

# -----------------------------
# INFO builder + het check + new helpers
# -----------------------------
def add_info_tags(info_field: str, tags: Dict[str, str]) -> str:
    if info_field in (".",""): info_field = ""
    items = [] if not info_field else info_field.split(";")
    for k,v in tags.items():
        if v is True or v == "1": items.append(k)
        else: items.append(f"{k}={v}")
    return ";".join([i for i in items if i]) or "."

def is_het(gt: Optional[str]) -> bool:
    if gt is None: return False
    g = gt.replace("|","/").split("/")
    return len(g)==2 and g[0] not in (".","") and g[1] not in (".","") and g[0] != g[1]

def child_has_alt(gt: Optional[str]) -> bool:
    """True if child GT has any ALT (>0). Tolerates './1','1/.','.|1','1|.'."""
    return any(a > 0 for a in (parse_gt_alleles(gt) or []))

def parent_has_alt(gt: Optional[str]) -> bool:
    """True if parent GT has any ALT (>0). Tolerates './1','1/.','.|1','1|.'."""
    return any(a > 0 for a in (parse_gt_alleles(gt) or []))

def haploid_missing_allele(gt: Optional[str]) -> bool:
    """True if GT is partially missing but includes an ALT (e.g., './1','1/.','.|1','1|.')."""
    if gt is None: return False
    parts = gt.replace("|","/").split("/")
    return (len(parts) == 2) and ("." in parts) and any(p not in ("0",".") for p in parts)

# -----------------------------
# FORMAT handling
# -----------------------------
def fetch_required_indices(fmt_keys: List[str], required: List[str]) -> Tuple[Dict[str,int], List[str]]:
    idx = {}; missing=[]
    for k in required:
        if k in fmt_keys: idx[k]=fmt_keys.index(k)
        else: missing.append(k)
    return idx, missing

def log_missing_required(missing: List[str], context: str, reported: Set[Tuple[str,...]], strict: bool, present_fmt: Optional[str]=None):
    key = tuple(sorted(missing))
    if key not in reported:
        reported.add(key)
        present = f" present={present_fmt}" if present_fmt else ""
        info(f"[ERROR] missing required FORMAT {missing} at {context}{present} "
             f"({'fatal' if strict else 'skipping'})")
    if strict:
        raise RuntimeError(f"Missing required FORMAT {missing} at {context} (present={present_fmt})")

# -----------------------------
# Region masking
# -----------------------------
class IntervalMask:
    def __init__(self): self.iv: Dict[str, List[Tuple[int,int]]] = {}
    def add(self, chrom: str, s: int, e: int): self.iv.setdefault(chrom, []).append((s,e))
    def build_from_bed_stream(self, fh: Iterable[str], name_hint: str=""):
        for line in fh:
            if not line.strip() or line.startswith("#"): continue
            p = line.split()
            if len(p)<3: continue
            chrom, s, e = p[0], p[1], p[2]
            try: s0, e0 = int(s), int(e)
            except ValueError: continue
            if s0>e0: s0,e0 = e0,s0
            self.add(chrom, s0+1, e0)
    def contains(self, chrom: str, pos: int) -> bool:
        for key in (chrom, chrom[3:] if chrom.startswith("chr") else "chr"+chrom):
            for s,e in self.iv.get(key, []):
                if s <= pos <= e: return True
        return False

def load_regions_mask(regions_path: str, verbose: bool=False) -> Optional[Tuple[IntervalMask, List[str]]]:
    if not regions_path: return None
    p = Path(regions_path); mask = IntervalMask(); loaded=[]
    def _load_stream(name: str, fh: Iterable[str]): mask.build_from_bed_stream(fh, name); loaded.append(name)
    if p.is_file() and str(p).endswith((".tar",".tar.gz",".tgz",".tar.bz2",".tbz2")):
        with tarfile.open(p, "r:*") as tf:
            for m in [m for m in tf.getmembers() if m.isreg()]:
                lower = m.name.lower()
                if lower.endswith(".bed") or lower.endswith(".bed.gz"):
                    try: data = tf.extractfile(m).read()
                    except Exception: continue
                    if lower.endswith(".gz"):
                        try: data = gzip.decompress(data)
                        except Exception: continue
                    fh = io.TextIOWrapper(io.BytesIO(data), encoding="utf-8")
                    _load_stream(Path(m.name).name.replace(".gz",""), fh)
    elif p.is_dir():
        for cand in p.rglob("*"):
            if not cand.is_file(): continue
            lower = cand.name.lower()
            if lower.endswith(".bed") or lower.endswith(".bed.gz"):
                try:
                    if lower.endswith(".gz"):
                        with gzip.open(cand,"rt") as fh: _load_stream(cand.name.replace(".gz",""), fh)
                    else:
                        with open(cand,"rt") as fh: _load_stream(cand.name, fh)
                except Exception: continue
    else:
        log(f"[regions] path not found/unsupported: {regions_path}", verbose); return None
    if not loaded: log("[regions] no BED files found; masking disabled.", verbose); return None
    if verbose: log(f"[regions] loaded {len(loaded)} mask file(s): {', '.join(loaded)}", True)
    return mask, loaded

# -----------------------------
# Denovo tagging (sex-aware) + scan with allele-level option
# -----------------------------
def tag_autosomal(f_gt: Optional[str], m_gt: Optional[str], p_gt: Optional[str]) -> bool:
    return (f_gt and m_gt and p_gt) and is_hom_ref(f_gt) and is_hom_ref(m_gt) and is_non_ref(p_gt)

def scan_trio_denovos(trio_gz: str, father_id: str, mother_id: str, child_id: str,
                      child_sex: Optional[str], par: ParIntervals, verbose: bool=False,
                      intersect_mode: str="allele", normalize_alleles_flag: bool=True
                      ) -> Tuple[List[str], Dict[str, Tuple[List[str], int, int, str, str]]]:
    """
    Return (header_lines, dict keyed by CHR:POS:REF:ALT -> (row, alt_index, norm_pos, norm_ref, norm_alt))
    If intersect_mode=='allele': split multi-allelic and require child carries that specific ALT.
    If 'locus': behave like original (no per-allele split), alt_index=-1.
    """
    denovo: Dict[str, Tuple[List[str], int, int, str, str]] = {}
    header_lines: List[str] = []
    reported_missing_sets: Set[Tuple[str,...]] = set()

    with gzip.open(trio_gz, "rt") as fh:
        name_to_idx = None
        for line in fh:
            if line.startswith("##"): header_lines.append(line); continue
            if line.startswith("#CHROM"):
                header_lines.append(line)
                cols = line.rstrip("\n").split("\t")
                samples = cols[9:]
                name_to_idx = {name: 9+i for i,name in enumerate(samples)}
                for sid in (father_id, mother_id, child_id):
                    if sid not in name_to_idx: raise ValueError(f"Sample {sid} not found in {trio_gz}")
                continue
            if name_to_idx is None: continue

            row = line.rstrip("\n").split("\t")
            chrom, pos_s, _vid, ref, alt_field = row[0], row[1], row[2], row[3], row[4]
            try: pos = int(pos_s)
            except ValueError: pos = 0

            # light screen
            if chrom.startswith("chrUn") or "_random" in chrom or "_alt" in chrom: continue
            # NOTE: do not early-screen homopolymer at multi-allelic level; we enforce per-allele later
            cbase = chrom_base(chrom)
            if child_sex is None and cbase in ("X","Y"): continue

            # need GT
            fmt = row[8]; fmt_keys = fmt.split(":")
            if "GT" not in fmt_keys:
                log_missing_required(["GT"], f"{trio_gz}:{chrom}:{pos_s}", reported_missing_sets, strict=False, present_fmt=":".join(fmt_keys))
                continue

            fi, mi, ci = name_to_idx[father_id], name_to_idx[mother_id], name_to_idx[child_id]
            f_gt, m_gt, p_gt = parse_gt(row[fi]), parse_gt(row[mi]), parse_gt(row[ci])

            # sex-aware denovo tagging (naive)
            is_dnv_locus = False
            if cbase in ("X","Y") and (child_sex or "").lower().startswith("m") and not par.in_par(chrom, pos):
                if cbase=="X": is_dnv_locus = (m_gt and p_gt and is_hom_ref(m_gt) and is_non_ref(p_gt))
                else:          is_dnv_locus = (f_gt and p_gt and is_hom_ref(f_gt) and is_non_ref(p_gt))
            else:
                is_dnv_locus = tag_autosomal(f_gt, m_gt, p_gt)
            if not is_dnv_locus: continue

            # per-allele split if requested
            alts = alt_field.split(",")
            if intersect_mode == "locus":
                npos, nref, nalt = (pos, ref, alts[0])
                if normalize_alleles_flag and not is_symbolic(nalt):
                    npos, nref, nalt = normalize_allele(pos, ref, alts[0])
                key = f"{chrom}:{npos}:{nref}:{nalt}"
                denovo[key] = (row, -1, npos, nref, nalt)
            else:
                child_alleles = parse_gt_alleles(p_gt) or []
                for i, alt in enumerate(alts, start=1):
                    if is_symbolic(alt): continue
                    if i not in child_alleles:  # child must carry this ALT
                        continue
                    npos, nref, nalt = (pos, ref, alt)
                    if normalize_alleles_flag:
                        npos, nref, nalt = normalize_allele(pos, ref, alt)
                    key = f"{chrom}:{npos}:{nref}:{nalt}"
                    denovo[key] = (row, i, npos, nref, nalt)

    log(f"[scan] {trio_gz} denovos: {len(denovo)}", verbose)
    return header_lines, denovo

# -----------------------------
# Metrics / sorting helpers
# -----------------------------
class Metrics:
    def __init__(self, child):
        self.child = child
        self.counts = Counter(scanned=0, denovo_tagged=0, region_kept=0, format_ok=0, sb_ok=0, final=0)
        self.dp = {"F":[], "M":[], "C":[]}; self.gq = {"F":[], "M":[], "C":[]}; self.ab=[]
        self.chrom_counts = Counter(); self.multi_allelic = 0; self.parent_leak = 0
    def bump(self, k): self.counts[k]+=1
    def add_dp_gq(self, fdp, mdp, kdp, fgq, mgq, kgq):
        if fdp is not None: self.dp["F"].append(fdp)
        if mdp is not None: self.dp["M"].append(mdp)
        if kdp is not None: self.dp["C"].append(kdp)
        if fgq is not None: self.gq["F"].append(fgq)
        if mgq is not None: self.gq["M"].append(mgq)
        if kgq is not None: self.gq["C"].append(kgq)
    def add_ab(self, val):
        if val is not None: self.ab.append(val)
    def as_dict(self):
        mids = lambda xs: statistics.median(xs) if xs else None
        total = sum(self.chrom_counts.values()) or 1
        pct = lambda n: round(100.0*n/total, 3)
        return {
            "child": self.child, "counts": dict(self.counts),
            "medians": {"DP_F": mids(self.dp["F"]), "DP_M": mids(self.dp["M"]), "DP_C": mids(self.dp["C"]),
                        "GQ_F": mids(self.gq["F"]), "GQ_M": mids(self.gq["M"]), "GQ_C": mids(self.gq["C"]),
                        "AB_C": mids(self.ab)},
            "pct": {"on_chrX": pct(self.chrom_counts.get("X",0)+self.chrom_counts.get("chrX",0)),
                    "on_chrY": pct(self.chrom_counts.get("Y",0)+self.chrom_counts.get("chrY",0)),
                    "on_PAR": pct(self.chrom_counts.get("PAR",0)),
                    "multi_allelic": pct(self.multi_allelic)},
            "parent_leak_sites": int(self.parent_leak),
        }

def header_contig_order(header_lines):
    order = {}; idx=0
    for h in header_lines:
        if h.startswith("##contig="):
            m = re.search(r"ID=([^,>]+)", h)
            if m: order[m.group(1)] = idx; idx+=1
    return order or None

def make_sort_key_fn(contig_order):
    def nat_key(chrom):
        parts = re.split(r'(\d+)', chrom)
        return [int(p) if p.isdigit() else p for p in parts]
    if not contig_order:
        return lambda chrom,pos: (nat_key(chrom), pos)
    return lambda chrom,pos: (contig_order.get(chrom, 10**9), pos)

# -----------------------------
# Trio readers & merge
# -----------------------------
def read_trios(family_file: str, default_child_sex: Optional[str]=None) -> List[Tuple[str,str,str,Optional[str]]]:
    trios=[]
    with open(family_file) as fh:
        for line in fh:
            line=line.strip()
            if not line or line.startswith("#"): continue
            parts = re.split(r"[,\s]+", line)
            if len(parts) < 3: continue
            head=[p.lower() for p in parts[:3]]
            if head[0].startswith("father") and head[1].startswith("mother"): continue
            sex=None
            if len(parts)>=4 and parts[3]: sex=parts[3].lower()
            elif default_child_sex: sex=default_child_sex.lower()
            trios.append((parts[0], parts[1], parts[2], sex))
    if not trios: raise RuntimeError("No trios found in family_file.")
    return trios

def read_trios_from_ped(ped_path: str, default_child_sex: Optional[str]=None) -> List[Tuple[str,str,str,Optional[str]]]:
    inds: Dict[str, Tuple[str,str,str]] = {}; present: Set[str] = set()
    with open(ped_path) as fh:
        for line in fh:
            line=line.strip()
            if not line or line.startswith("#"): continue
            p=line.split()
            if len(p)<6: continue
            _fid, iid, pid, mid, sex, _pheno = p[:6]
            present.add(iid); inds[iid]=(pid,mid,sex)
    trios=[]
    for iid,(pid,mid,sex) in inds.items():
        if not pid or pid=="0" or not mid or mid=="0": continue
        if pid not in present or mid not in present: continue
        sex_str=None
        if sex=="1": sex_str="male"
        elif sex=="2": sex_str="female"
        elif default_child_sex: sex_str=default_child_sex.lower()
        trios.append((pid,mid,iid,sex_str))
    if not trios: raise RuntimeError("No complete trios found in PED.")
    return trios

def merge_trios(family_trios, ped_trios, prefer="ped"):
    by_child: Dict[str, Tuple[str,str,str,Optional[str],str]] = {}
    def put(trio, source):
        f,m,c,sex = trio
        if c not in by_child: by_child[c]=(f,m,c,sex,source); return
        f0,m0,_,sex0,src0 = by_child[c]
        if (f0!=f) or (m0!=m):
            info(f"[merge] conflicting parents for child {c}: {src0}({f0},{m0}) vs {source}({f},{m}) -> pref {prefer}")
            if source==prefer: by_child[c]=(f,m,c, sex if sex else sex0, source)
        else:
            if (sex0 is None) and sex: by_child[c]=(f0,m0,c,sex,src0)
            elif sex0 and sex and sex0!=sex:
                info(f"[merge] conflicting sex for child {c}: {src0}({sex0}) vs {source}({sex}) -> pref {prefer}")
                if source==prefer: by_child[c]=(f0,m0,c,sex,source)
    if prefer=="ped":
        for t in ped_trios: put(t,"ped")
        for t in family_trios: put(t,"family")
    else:
        for t in family_trios: put(t,"family")
        for t in ped_trios: put(t,"ped")
    return [(f,m,c,sex) for (f,m,c,sex,_) in by_child.values()]

# -----------------------------
# Main trio processing
# -----------------------------
def process_trio(args, f, m, c, child_sex):
    vbool = args.verbose > 0
    par = ParIntervals(args.par_bed, verbose=vbool)

    effective_sex = (child_sex or args.default_child_sex)
    if effective_sex: info(f"[sex] child {c}: {'provided' if child_sex else 'defaulted'} {effective_sex}")
    else: info(f"[sex] child {c}: sex unknown → autosomes-only (chr1–22)")

    c1_base = Path(args.caller1_vcf)
    c2_base = Path(args.caller2_vcf) if args.caller2_vcf else None
    c1 = c1_base / f"{c}.trio.caller1.vcf.gz" if c1_base.is_dir() else c1_base
    c2 = c2_base / f"{c}.trio.caller2.vcf.gz" if (c2_base and c2_base.is_dir()) else (c2_base if c2_base else None)

    if not c1.exists(): raise FileNotFoundError(f"Missing required caller1 VCF: {c1}")
    dual_mode=False
    if c2_base:
        if c2 and c2.exists():
            dual_mode=True; info(f"[mode] two-caller intersection (caller1 ∩ caller2) for child {c}")
        else:
            info(f"[mode] caller2 not found ({c2}); caller1-only for child {c}")
    if not c2_base: info(f"[mode] caller2 not provided; caller1-only for child {c}")

    if args.dry_run:
        info(f"[dry-run] would process child={c} with c1={c1} c2={c2} out={args.out_dir}")
        print(str(Path(args.out_dir)/f"{c}.final.vcf.gz")); return

    # Scan caller(s)
    c1_hdr, c1_dnvs = scan_trio_denovos(str(c1), f, m, c, effective_sex, par, vbool,
                                        intersect_mode=args.intersect_mode, normalize_alleles_flag=args.normalize_alleles)
    if dual_mode:
        c2_hdr, c2_dnvs = scan_trio_denovos(str(c2), f, m, c, effective_sex, par, vbool,
                                            intersect_mode=args.intersect_mode, normalize_alleles_flag=args.normalize_alleles)
        inter_keys = set(c1_dnvs.keys()) & set(c2_dnvs.keys())
        log(f"[intersect] shared keys: {len(inter_keys)}", vbool)
    else:
        inter_keys = set(c1_dnvs.keys())

    header_lines = c1_hdr if c1_hdr else (c2_hdr if dual_mode else [])
    name_to_idx = {}
    for l in header_lines:
        if l.startswith("#CHROM"):
            cols = l.rstrip("\n").split("\t")
            name_to_idx = {name: 9+i for i,name in enumerate(cols[9:])}
            break
    for sid in (f,m,c):
        if sid not in name_to_idx: raise RuntimeError(f"Sample {sid} not present in header.")

    mask_result = load_regions_mask(args.regions, verbose=vbool) if args.regions else None
    if mask_result:
        region_mask, loaded_list = mask_result; info(f"[regions] final output WILL be region-filtered ({', '.join(loaded_list)})")
    else:
        region_mask=None
        info("[regions] not requested → NOT region-filtered" if not args.regions else "[regions] requested, but no masks found → NOT region-filtered")

    required_for_filters = [k for k in (args.require_format_keys or "").split(",") if k] or ["GT","DP","AD","GQ"]
    strict_mode = bool(args.strict_format)
    reported_missing_sets: Set[Tuple[str,...]] = set()

    # Output paths
    out_dir = Path(args.out_dir); out_dir.mkdir(parents=True, exist_ok=True)
    final_gz = out_dir / f"{c}.final.vcf.gz"
    tmp_vcf = out_dir / f"{c}.final.tmp.vcf"  # plain text; will bgzip if pysam available
    sites_tsv_path = out_dir / f"{c}.sites.tsv"
    sites_tsv_fh = open(sites_tsv_path, "wt") if args.emit_sites_tsv else None
    if sites_tsv_fh:
        sites_tsv_fh.write("\t".join(["CHROM","POS","REF","ALT","F_GT","M_GT","C_GT","F_DP","M_DP","C_DP","F_GQ","M_GQ","C_GQ","C_AB"]) + "\n")

    metrics = Metrics(c)
    contig_order = header_contig_order(header_lines)
    sort_key_fn = make_sort_key_fn(contig_order)

    # --- write plain VCF first ---
    with open(tmp_vcf, "wt") as out:
        version_line = "##HAT-FLEX_version=1.0\n"
        wrote_version=False

        # Detect whether key INFO headers already exist
        have_denovo = any(h.startswith("##INFO=<ID=DENOVO") for h in header_lines)
        have_sexctx = any(h.startswith("##INFO=<ID=SEXCTX") for h in header_lines)
        have_irrel  = any(h.startswith("##INFO=<ID=IRRELPARENT") for h in header_lines)

        for h in header_lines:
            if h.startswith("##"):
                out.write(h)
                if h.startswith("##fileformat"): out.write(version_line); wrote_version=True
        if not wrote_version: out.write(version_line)

        # Add missing INFO declarations used by this pipeline
        if not have_denovo:
            out.write('##INFO=<ID=DENOVO,Number=0,Type=Flag,Description="De novo candidate at trio locus">\n')
        if not have_sexctx:
            out.write('##INFO=<ID=SEXCTX,Number=1,Type=String,Description="Genomic sex context: Autosomal, PAR, MaleX, FemaleX, MaleY">\n')
        if not have_irrel:
            out.write('##INFO=<ID=IRRELPARENT,Number=1,Type=String,Description="Parent sample ignored for haploid locus: F or M">\n')
        # Unusual haploid signals
        out.write('##INFO=<ID=HAPUNUSUAL,Number=0,Type=Flag,Description="Unusual pattern at haploid locus (male non-PAR X/Y): het-like GT or mid-range AB">\n')
        out.write('##FILTER=<ID=HaploidUnusual,Description="Child haploid locus shows het-like GT or suspicious mid-range allele balance">\n')

        # Simple cluster annotation (PASS-only; present only on PASS records)
        out.write('##INFO=<ID=CLUSTER,Number=1,Type=String,Description="PASS-only; yes/no if within --cluster-window-bp among PASS de novos">\n')

        for h in header_lines:
            if h.startswith("#CHROM"):
                cols = h.rstrip("\n").split("\t")
                trio_header = cols[:9] + [f, m, c]
                out.write("\t".join(trio_header) + "\n"); break

        records: List[Tuple[Tuple, str, List[str]]] = []
        # Per-chromosome sliding window for clustering (PASS-only).
        # We store (position, record_index) in the deque so we can retro-update previous PASS records.
        cluster_state  = defaultdict(lambda: deque())

        def _bump_many_to_cluster_yes(records_list, idx_list):
            """Retroactively set CLUSTER=yes on previously stored PASS records."""
            for prev_idx in idx_list or []:
                if prev_idx is None or prev_idx < 0 or prev_idx >= len(records_list):
                    continue
                key2, line, reasons = records_list[prev_idx]
                parts = line.split("\t")
                if len(parts) < 8:
                    continue
                info_field = parts[7]
                if "CLUSTER=yes" in info_field:
                    pass
                elif "CLUSTER=no" in info_field:
                    info_field = info_field.replace("CLUSTER=no", "CLUSTER=yes")
                else:
                    info_field = (info_field if info_field != "." else "") + (";" if info_field and info_field != "." else "") + "CLUSTER=yes"
                parts[7] = info_field or "."
                records_list[prev_idx] = (key2, "\t".join(parts), reasons)

        # --------- sort de novos by (chrom, npos) ---------
        def _get_dnvo(k):
            return c1_dnvs.get(k) or (c2_dnvs.get(k) if dual_mode else None)
        ordered_keys = sorted(
            inter_keys,
            key=lambda k: sort_key_fn(_get_dnvo(k)[0][0], _get_dnvo(k)[2])  # (chrom, npos)
        )
        # ---------------------------------------------------

        for key in ordered_keys:
            metrics.bump("denovo_tagged")
            r1 = _get_dnvo(key)
            if r1 is None: continue
            row, alt_index, npos, nref, nalt = r1
            chrom = row[0]; pos_s = row[1]; ref = row[3]; alt_field = row[4]
            try: orig_pos = int(pos_s)
            except ValueError: orig_pos = npos

            filt_reasons: List[str] = []
            metrics.bump("scanned")

            # count multi-allelic loci
            if "," in alt_field:
                metrics.multi_allelic += 1

            if is_symbolic(nalt): filt_reasons.append("Symbolic")

            if (effective_sex is None) and chrom_base(chrom) in ("X","Y"):
                filt_reasons.append("SexSkip")

            # region mask on normalized pos
            region_kept_here = True
            if region_mask and region_mask.contains(chrom, npos):
                filt_reasons.append("Masked"); region_kept_here=False
            if region_kept_here: metrics.bump("region_kept")

            if is_bad_homopolymer_param(nref, nalt, args.max_homopolymer):
                filt_reasons.append("Homopolymer")

            fmt = row[8]; fmt_keys = fmt.split(":")
            idx_map, missing = fetch_required_indices(fmt_keys, required_for_filters)
            if missing:
                log_missing_required(missing, f"{chrom}:{npos}", reported_missing_sets, strict_mode, present_fmt=":".join(fmt_keys))
                filt_reasons.append("MissingFORMAT")
            else:
                metrics.bump("format_ok")

            # Trio fields
            try:
                s_f = row[name_to_idx[f]]; s_m = row[name_to_idx[m]]; s_c = row[name_to_idx[c]]
            except Exception:
                filt_reasons.append("SampleCols"); s_f=s_m=s_c=None

            # SB child
            sb_child = parse_format_sb(s_c or "", fmt_keys) if s_c else None
            if not alt_on_both_strands(sb_child, args.sb_min_each): filt_reasons.append("SB")
            if (args.sb_fisher_pmax is not None) and sb_child and len(sb_child)==4:
                a,b,cF,d = sb_child
                try:
                    p = fisher_right_tail(a,b,cF,d)  # right-tail by design
                    if p > args.sb_fisher_pmax: filt_reasons.append("SBP")
                except Exception: pass

            # child genotype must include this allele if allele-level mode
            child_gt = parse_gt(s_c or "")
            child_alleles = parse_gt_alleles(child_gt) or []
            if args.intersect_mode=="allele" and alt_index>0 and (alt_index not in child_alleles):
                filt_reasons.append("ChildNoAllele")

            # exact haploid context at this locus
            child_hap = child_is_haploid(effective_sex, chrom, npos, par)

            # ---------------------------
            # Haploid-aware thresholds & relevant-parent selection
            # ---------------------------
            eff_child_dp = args.haploid_depth_value if child_hap else args.depth_value
            eff_fd = args.depth_value
            eff_md = args.depth_value

            req_f, req_m = True, True
            haploid_context = child_hap  # precise
            if haploid_context:
                cbase = chrom_base(chrom)
                if cbase == "X":
                    req_f, req_m = False, True
                    eff_md = args.haploid_depth_value
                elif cbase == "Y":
                    req_f, req_m = True, False
                    eff_fd = args.haploid_depth_value

            # Numeric fields
            iDP = idx_map.get("DP"); iAD=idx_map.get("AD"); iGQ=idx_map.get("GQ")
            fd = md = kd = fg = mg = kg = None
            fa = ma = ka = 0
            try:
                fS = s_f.split(":") if s_f else []
                mS = s_m.split(":") if s_m else []
                kS = s_c.split(":") if s_c else []
                fd = int(fS[iDP]) if (iDP is not None and iDP < len(fS)) else None
                md = int(mS[iDP]) if (iDP is not None and iDP < len(mS)) else None
                kd = int(kS[iDP]) if (iDP is not None and iDP < len(kS)) else None
                def alt_count(val, ai):
                    if not val or "," not in val: return 0
                    parts = val.split(",")
                    if ai <= 0:  # fall back to first ALT
                        return int(parts[1]) if len(parts)>1 else 0
                    return int(parts[ai]) if len(parts)>ai else 0
                if iAD is not None:
                    fa = alt_count(fS[iAD] if iAD < len(fS) else "", alt_index)
                    ma = alt_count(mS[iAD] if iAD < len(mS) else "", alt_index)
                    ka = alt_count(kS[iAD] if iAD < len(kS) else "", alt_index)
                fg = int(fS[iGQ]) if (iGQ is not None and iGQ < len(fS)) else None
                mg = int(mS[iGQ]) if (iGQ is not None and iGQ < len(mS)) else None
                kg = int(kS[iGQ]) if (iGQ is not None and iGQ < len(kS)) else None
            except Exception:
                filt_reasons.append("BadFORMAT")

            # Child presence check (GT-tolerant for haploid)
            if child_hap:
                if not child_has_alt(child_gt):
                    filt_reasons.append("ChildRef")
            else:
                if not is_non_ref(child_gt):
                    filt_reasons.append("ChildRef")

            # ----- Parent ALT-leak (haploid override + GT-based option) -----
            haploid_max = args.haploid_max_parent_alt_reads if args.haploid_max_parent_alt_reads is not None else args.max_parent_alt_reads
            eff_max_f = args.max_parent_alt_reads
            eff_max_m = args.max_parent_alt_reads
            if haploid_context:
                cbase = chrom_base(chrom)
                if cbase == "X":
                    eff_max_m = haploid_max  # mother relevant on X
                elif cbase == "Y":
                    eff_max_f = haploid_max  # father relevant on Y

            # (A) AD-based leak checks
            if iAD is not None:
                if (req_f and fa > eff_max_f) or (req_m and ma > eff_max_m):
                    filt_reasons.append("ParentAlt")
            if iAD is not None and ka < args.min_child_alt_reads:
                filt_reasons.append("ChildAlt")

            # (B) GT-based leak checks for relevant parent(s) if enabled
            if args.parent_gt_alt_triggers_parentalt:
                if req_f and parent_has_alt(parse_gt(s_f or "")):
                    filt_reasons.append("ParentAlt")
                if req_m and parent_has_alt(parse_gt(s_m or "")):
                    filt_reasons.append("ParentAlt")

            # ----- DP/GQ checks (only relevant parent(s); child uses haploid-aware DP) -----
            if iDP is not None and iGQ is not None:
                needed_vals = []
                if req_f: needed_vals += [fd, fg]
                if req_m: needed_vals += [md, mg]
                needed_vals += [kd, kg]

                if any(v is None for v in needed_vals):
                    filt_reasons.append("MissingDPGQ")
                else:
                    if req_f:
                        if fd < eff_fd: filt_reasons.append("LowDP")
                        if fg <= args.gq_value: filt_reasons.append("LowGQ")
                        if iAD is not None and fa > eff_max_f: filt_reasons.append("ParentAlt")
                    if req_m:
                        if md < eff_md: filt_reasons.append("LowDP")
                        if mg <= args.gq_value: filt_reasons.append("LowGQ")
                        if iAD is not None and ma > eff_max_m: filt_reasons.append("ParentAlt")
                    if kd < eff_child_dp: filt_reasons.append("LowDP")
                    if kg <= args.gq_value: filt_reasons.append("LowGQ")

            # AB checks (global vs haploid-specific floor)
            ab_val = (float(ka)/kd) if (iAD is not None and kd and kd>0) else None
            if iDP is not None and iAD is not None and kd and kd>0:
                if child_hap:
                    if float(ka)/kd < args.haploid_ab_min:
                        filt_reasons.append("LowAB_haploid")
                else:
                    if float(ka)/kd < args.ab_min:
                        filt_reasons.append("LowAB")

            # --- Flag unusual haploid pattern (MaleX/MaleY non-PAR): het-like or mid-range AB ---
            haploid_unusual = False
            if child_hap:
                missing_one = haploid_missing_allele(child_gt)
                if is_het(child_gt):
                    haploid_unusual = True
                elif ab_val is not None and args.haploid_suspicious_ab_low <= ab_val <= args.haploid_suspicious_ab_high:
                    if not (args.haploid_suspicious_skip_missing_allele and missing_one):
                        haploid_unusual = True
            if haploid_unusual:
                filt_reasons.append("HaploidUnusual")

            metrics.add_ab(ab_val); metrics.add_dp_gq(fd, md, kd, fg, mg, kg)
            if par.in_par(chrom, npos): metrics.chrom_counts["PAR"] += 1
            metrics.chrom_counts[chrom] += 1

            # Build output row with normalized POS/REF/ALT, but keep original INFO/FORMAT
            out_cols = list(row)
            out_cols[1] = str(npos); out_cols[3] = nref; out_cols[4] = nalt

            # --- SEXCTX with FemaleX ---
            sex_l = (effective_sex or "").lower()
            in_par = par.in_par(chrom, npos)
            cbase = chrom_base(chrom)
            if in_par:
                sexctx = "PAR"
            elif cbase == "X" and sex_l.startswith("f"):
                sexctx = "FemaleX"
            elif cbase == "X" and child_is_haploid(effective_sex, chrom, npos, par):
                sexctx = "MaleX"
            elif cbase == "Y" and child_is_haploid(effective_sex, chrom, npos, par):
                sexctx = "MaleY"
            else:
                sexctx = "Autosomal"

            # Start INFO tags (clustering will be added ONLY for PASS sites below)
            tags = {"DENOVO": "1", "SEXCTX": sexctx}

            if ("SB" not in filt_reasons) and ("SBP" not in filt_reasons): metrics.bump("sb_ok")

            passes = (len(filt_reasons) == 0)

            # --- Compute clustering ONLY for PASS sites (and update deque with PASS positions only) ---
            is_cluster = False
            idxs_to_upgrade = []
            if passes:
                # Single window
                win = int(args.cluster_window_bp)
                dq = cluster_state[chrom]
                # Drop positions older than window
                while dq and (npos - dq[0][0]) > win:
                    dq.popleft()
                # This site is in a cluster if there is at least (cluster_min_count-1) PASS site(s) already in window
                is_cluster = len(dq) >= max(0, int(args.cluster_min_count) - 1)
                # If clustered, remember to flip ALL prior PASS in-window to CLUSTER=yes as well.
                if is_cluster and len(dq) > 0:
                    idxs_to_upgrade = [rec_idx for (_pos, rec_idx) in dq]
                tags.update({"CLUSTER": "yes" if is_cluster else "no"})

            info_field = add_info_tags(out_cols[7], tags)

            # Annotate which parent was ignored (if exactly one is ignored)
            if (req_f ^ req_m):
                irrelevant = "M" if req_f and not req_m else "F"
                info_field = add_info_tags(info_field, {"IRRELPARENT": irrelevant})
            # Mark HAPUNUSUAL in INFO if flagged
            if haploid_unusual:
                info_field = add_info_tags(info_field, {"HAPUNUSUAL":"1"})
            out_cols[7] = info_field

            # Hide irrelevant parent's FORMAT field in the written VCF
            sf_out = (s_f if req_f else ".")
            sm_out = (s_m if req_m else ".")
            trio_fixed = out_cols[:9] + [sf_out or ".", sm_out or ".", s_c or "."]

            if passes:
                metrics.bump("final")
                if iAD is not None and ((fa>0 and req_f) or (ma>0 and req_m)): metrics.parent_leak += 1

            if args.emit_sites_tsv:
                FGT, MGT, CGT = parse_gt(s_f or "") or ".", parse_gt(s_m or "") or ".", parse_gt(s_c or "") or "."
                FDP_raw = str(fd) if fd is not None else "."
                MDP_raw = str(md) if md is not None else "."
                KDP = str(kd) if kd is not None else "."
                FGQ_raw = str(fg) if fg is not None else "."
                MGQ_raw = str(mg) if mg is not None else "."
                KGQ = str(kg) if kg is not None else "."
                CAB = f"{ab_val:.3f}" if ab_val is not None else "."

                # Mask irrelevant parent's fields in TSV output
                FGT = FGT if req_f else "."
                MGT = MGT if req_m else "."
                FDP = FDP_raw if req_f else "."
                MDP = MDP_raw if req_m else "."
                FGQ = FGQ_raw if req_f else "."
                MGQ = MGQ_raw if req_m else "."

                sites_tsv_fh.write("\t".join([chrom, str(npos), nref, nalt, FGT, MGT, CGT, FDP, MDP, KDP, FGQ, MGQ, KGQ, CAB]) + "\n")

            if args.keep_failures or passes:
                filt_val = "PASS" if passes else ";".join(sorted(set(filt_reasons)))
                trio_fixed[6] = filt_val
                key2 = sort_key_fn(chrom, npos)
                # Append the record line now
                records.append((key2, "\t".join(trio_fixed), filt_reasons))

                # If this PASS record formed a cluster, retroactively mark ALL prior PASS in-window.
                if passes and is_cluster and idxs_to_upgrade:
                    _bump_many_to_cluster_yes(records, idxs_to_upgrade)

                # Finally, if PASS, add this site to the window with its index for future retro-updates.
                if passes:
                    this_idx = len(records) - 1
                    dq = cluster_state[chrom]
                    dq.append((npos, this_idx))

            if args.verbose >= 2 and (len(records) % 10000 == 0) and len(records)>0:
                info(f"[progress] child={c} records_written={len(records)}")

        records.sort(key=lambda x: x[0])
        for _, line, _reasons in records:
            out.write(line + "\n")

    if sites_tsv_fh: sites_tsv_fh.close()

    # --- compress & index ---
    if HAVE_PYSAM:
        try:
            pysam.tabix_compress(str(tmp_vcf), str(final_gz), force=True)
            if args.write_index != "none":
                pysam.tabix_index(str(final_gz), preset="vcf", force=True, csi=(args.write_index=="csi"))
        except Exception as e:
            info(f"[index] pysam bgzip/index failed ({e}); falling back to gzip (no index).")
            with open(tmp_vcf,"rt") as inp, gzip.open(final_gz,"wt") as outgz:
                shutil.copyfileobj(inp, outgz)
    else:
        if args.write_index != "none":
            info("[index] --write-index requested but pysam not available; writing gzip without index.")
        with open(tmp_vcf,"rt") as inp, gzip.open(final_gz,"wt") as outgz:
            shutil.copyfileobj(inp, outgz)

    try: os.remove(tmp_vcf)
    except Exception: pass

    # Metrics file
    if args.metrics:
        m = metrics.as_dict()
        outp = Path(args.metrics)
        if outp.suffix.lower()==".json":
            with open(outp, "at" if outp.exists() else "wt") as fh: fh.write(json.dumps(m)+"\n")
        else:
            write_header = (not outp.exists())
            with open(outp, "at") as fh:
                if write_header:
                    fh.write("\t".join(["child","scanned","denovo_tagged","region_kept","format_ok","sb_ok","final",
                                        "med_DP_F","med_DP_M","med_DP_C","med_GQ_F","med_GQ_M","med_GQ_C","med_AB_C",
                                        "pct_chrX","pct_chrY","pct_PAR","pct_multi","parent_leak_sites"]) + "\n")
                fh.write("\t".join(map(str, [
                    m["child"], m["counts"]["scanned"], m["counts"]["denovo_tagged"], m["counts"]["region_kept"],
                    m["counts"]["format_ok"], m["counts"]["sb_ok"], m["counts"]["final"],
                    m["medians"]["DP_F"], m["medians"]["DP_M"], m["medians"]["DP_C"],
                    m["medians"]["GQ_F"], m["medians"]["GQ_M"], m["medians"]["GQ_C"], m["medians"]["AB_C"],
                    m["pct"]["on_chrX"], m["pct"]["on_chrY"], m["pct"]["on_PAR"], m["pct"]["multi_allelic"],
                    m["parent_leak_sites"]
                ])) + "\n")

    # Run manifest
    manifest = {
        "timestamp": time.strftime("%Y-%m-%dT%H:%M:%S%z"),
        "args": {k:(str(v) if isinstance(v, Path) else v) for k,v in vars(args).items()},
        "inputs": {"caller1": str(c1), "caller2": (str(c2) if c2 else None),
                   "family_file": args.family_file, "ped": args.ped, "regions": args.regions, "par_bed": args.par_bed},
        "outputs": {"final_vcf": str(final_gz),
                    "sites_tsv": (str(sites_tsv_path) if args.emit_sites_tsv else None)}
    }
    with open(Path(args.out_dir)/f"{c}.run_manifest.json","wt") as fh:
        json.dump(manifest, fh, indent=2)

    info(f"[final] wrote {final_gz}")
    print(str(final_gz))

# -----------------------------
# CLI
# -----------------------------
def main():
    epilog = r"""
Notes
-----
* If child sex is unknown, chrX/chrY are skipped for that child.
* Male non-PAR X/Y: haploid checks use --haploid-depth-value (=5 default) for child and relevant parent.
* Parent ALT leak can be detected by AD and (optionally) by GT (./1, 1/., .|1, 1|.) via
  --parent-gt-alt-triggers-parentalt / --no-parent-gt-alt-triggers-parentalt (applies only to relevant parent on haploid loci).
* Haploid "unusual" mid-range AB can be skipped when GT is './ALT' or 'ALT/.' at haploid loci via
  --haploid-suspicious-skip-missing-allele / --no-haploid-suspicious-skip-missing-allele.
* CLUSTER is a PASS-only annotation: only PASS records receive CLUSTER=yes|no based on nearby PASS de novos.
"""

    from argparse import ArgumentDefaultsHelpFormatter

    ap = argparse.ArgumentParser(
        description="Welcome to HAT-FLEX for calling de novo variants. Check out TNTurnerLab GitHub to learn more.",
        formatter_class=ArgumentDefaultsHelpFormatter,
        epilog=epilog
    )

    # Trio sources
    ap.add_argument("--family-file", default=None, help="This file contains the fatherName,motherName,childName and is optional. If not used the pedigree file will be used instead.")
    ap.add_argument("--ped", default=None, help="Standard pedigree file of fid, iid, pid, mid, sex, pheno")
    ap.add_argument("--prefer-source", choices=["ped","family"], default="ped", help="If both exist, which one is preferred by the user the family-file or the pedigree file")

    # Callers and intersection behavior
    ap.add_argument("--caller1-vcf", required=True, help="Path to Caller 1 VCF")
    ap.add_argument("--caller2-vcf", default=None, help="Path to Caller 2 VCF")
    ap.add_argument("--intersect-mode", choices=["locus","allele"], default="allele", help="Type of intersection; exact=allele or position=locus")
    ap.add_argument("--normalize-alleles", action="store_true", help="option to normalize alleles")

    # Regions / PAR
    ap.add_argument("--regions", default=None, help="Directory of regions (as compressed or uncompressed bed files) you want masked in the output (e.g., segmental duplications)")
    ap.add_argument("--par-bed", default=None, help="bed file containing coordinates of pseudoautosomal regions on the X and Y chromosome. If this is not provided, the default is b38 of the human genome.")

    # Genotype/quality filters
    ap.add_argument("--gq-value", type=int, default=20, help="Minimum GQ desired for genotypes. Default is 20.")
    ap.add_argument("--depth-value", type=int, default=10, help="Minimum DP desired for site. Default is 10.")
    ap.add_argument("--haploid-depth-value", type=int, default=5,
        help="DP minimum when the child is haploid (male non-PAR X/Y). Also applied to the relevant parent at haploid loci. Default is 5.")
    ap.add_argument("--ab-min", type=float, default=0.25, help="Minimum AB desired for alternate allele. Default is 0.25")
    ap.add_argument("--haploid-ab-min", type=float, default=0.85, help="Minimum AB desired for alternate allele in the haploid state. Default is 0.85")
    ap.add_argument("--haploid-suspicious-ab-low", type=float, default=0.35, help="Suspicious low AB for the haploid state. Default is 0.35")
    ap.add_argument("--haploid-suspicious-ab-high", type=float, default=0.65, help="Suspicious high AB for the haploid state. Default is 0.65")
    ap.add_argument("--default-child-sex", choices=["male","female"], default=None, help="Child Sex. If not provided, only the autosomes will be run")
    ap.add_argument("--require-format-keys", default="GT,DP,AD,GQ", help="Format keys required for consideration in the filtering scheme.")
    ap.add_argument("--strict-format", action="store_true", help="If this is provided than all the format keys in --require-format-keys muts be present in the variant line for it to be considered.")

    # Strand bias / read-count knobs
    ap.add_argument("--sb-min-each", type=int, default=1, help="SB minimum optional.")
    ap.add_argument("--sb-fisher-pmax", type=float, default=None,
        help="Right-tail Fisher exact p-value maximum for strand bias on ALT reads (set to, e.g., 0.001).")
    ap.add_argument("--min-child-alt-reads", type=int, default=0, help="minimum number of reads the alternate allele must be seen in the child.")
    ap.add_argument("--max-parent-alt-reads", type=int, default=0,
        help="Maximum alternate allele reads allowed in each parent (AD-based). Overridden for the relevant parent at haploid loci by --haploid-max-parent-alt-reads if provided.")
    ap.add_argument("--haploid-max-parent-alt-reads", type=int, default=None,
        help="If set, overrides --max-parent-alt-reads for the relevant parent at haploid loci (MaleX→mother; MaleY→father).")
    ap.add_argument("--max-homopolymer", type=int, default=10, help="Maximum homopolymer length of A or T. Default is 10.")

    # NEW toggles (with on/off pairs and defaults ON)
    ap.add_argument("--parent-gt-alt-triggers-parentalt", dest="parent_gt_alt_triggers_parentalt", action="store_true",
                    help="If set, relevant parent's GT containing alternate allele triggers ParentAlt even if AD missing.")
    ap.add_argument("--no-parent-gt-alt-triggers-parentalt", dest="parent_gt_alt_triggers_parentalt", action="store_false",
                    help="Disable GT-based ParentAlt for the relevant parent.")
    ap.set_defaults(parent_gt_alt_triggers_parentalt=True)

    ap.add_argument("--haploid-suspicious-skip-missing-allele", dest="haploid_suspicious_skip_missing_allele", action="store_true",
                    help="Skip HaploidUnusual when child haploid GT is './ALT' or 'ALT/.' (partial-missing) at haploid loci.")
    ap.add_argument("--no-haploid-suspicious-skip-missing-allele", dest="haploid_suspicious_skip_missing_allele", action="store_false",
                    help="Do not skip HaploidUnusual for partially missing haploid GTs.")
    ap.set_defaults(haploid_suspicious_skip_missing_allele=True)

    # Clustering (INFO-only; simple PASS-based window)
    ap.add_argument("--cluster-window-bp", type=int, default=100,
        help="Clustering window (bp) for PASS variants; sets INFO/CLUSTER on PASS records to yes/no based on window content.")
    ap.add_argument("--cluster-min-count", type=int, default=2,
        help="Minimum PASS de novos within the window (including this site) to call CLUSTER=yes (default 2 = requires ≥1 other PASS site).")

    # Outputs / audit / indexing
    ap.add_argument("--out-dir", default="/path/to/outdir", help="output directory for the results.")
    ap.add_argument("--metrics", default=None, help="optional metrics file to store critical metrics.")
    ap.add_argument("--emit-sites-tsv", action="store_true", help="optional tsv file of sites")
    ap.add_argument("--keep-failures", action="store_true", help="option to keep the failing sites in the VCF with the listing of why they failed to be a true de novo")
    ap.add_argument("--write-index", choices=["none","tbi","csi"], default="none", help="option to write a tabix-index for the output vcf")

    ap.add_argument("--dry-run", action="store_true", help="dry run of the analysis")
    ap.add_argument("-v", "--verbose", action="count", default=0, help="verbose messaging of processing steps")

    args = argparse.Namespace(**vars(ap.parse_args()))
    Path(args.out_dir).mkdir(parents=True, exist_ok=True)

    if not Path(args.caller1_vcf).exists():
        print(f"ERROR: {args.caller1_vcf} not found", file=sys.stderr); sys.exit(1)
    if args.caller2_vcf and not Path(args.caller2_vcf).exists():
        info(f"[mode] caller2 base path not found: {args.caller2_vcf} → caller1-only mode will be used")

    _ = ParIntervals(args.par_bed, verbose=(args.verbose>0))

    # Gather trios
    family_trios=[]; ped_trios=[]
    if args.family_file:
        if Path(args.family_file).exists():
            family_trios = read_trios(args.family_file, args.default_child_sex)
            info(f"[trios] loaded {len(family_trios)} from family file")
        elif not args.ped:
            print(f"ERROR: --family-file {args.family_file} not found and no --ped provided.", file=sys.stderr); sys.exit(1)
        else:
            info(f"[trios] family file not found: {args.family_file} -> relying on PED")
    if args.ped:
        if Path(args.ped).exists():
            ped_trios = read_trios_from_ped(args.ped, args.default_child_sex)
            info(f"[trios] loaded {len(ped_trios)} from PED")
        elif not family_trios:
            print(f"ERROR: --ped {args.ped} not found and no valid --family-file provided.", file=sys.stderr); sys.exit(1)
        else:
            info(f"[trios] PED not found: {args.ped} -> relying on family file only")
    if not family_trios and not ped_trios:
        print("ERROR: provide at least one of --family-file or --ped (and ensure it exists).", file=sys.stderr); sys.exit(1)

    trios = merge_trios(family_trios, ped_trios, prefer=args.prefer_source)
    info(f"[trios] merged total: {len(trios)}")

    if args.dry_run:
        info(f"[dry-run] caller1={args.caller1_vcf} caller2={args.caller2_vcf} out={args.out_dir} regions={args.regions} pysam={'yes' if HAVE_PYSAM else 'no'}")
        for father, mother, child, sex in trios:
            info(f"[dry-run] trio F={father} M={mother} C={child} sex={sex or 'unknown'} mode={args.intersect_mode} normalize={args.normalize_alleles}")
            print(str(Path(args.out_dir)/f"{child}.final.vcf.gz"))
        return

    for father, mother, child, sex in trios:
        process_trio(args, father, mother, child, sex)

    log("[DONE] all trios processed.", args.verbose>0)

if __name__ == "__main__":
    main()

