#!/usr/bin/env python
"""Translate Panaroo reference genes and merge annotation search results.

  python merge_annotations.py translate pan_genome_reference.fa proteins.faa translation_qc.tsv
  python merge_annotations.py merge --rtab gene_presence_absence2.Rtab --meta isolate_inhibi_category.txt \
         --kegg full_kegg.txt --swissprot best_hits.tsv --out pangenome_annotation.tsv

Expects res/{uniref90.tsv,cog.tsv,pfam.domtbl,kofam_detail.tsv} and db/{cog-20.cog.csv,
cog-20.def.tab,fun-20.tab,pfam2go,ko_list} from annotate_pangenome.sh.
"""
import argparse, io, re, sys, time
import numpy as np, pandas as pd, requests

# ---------------------------------------------------------------- translate
def translate(fa, faa, qc):
    b = "TCAG"; aa = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
    code = {x + y + z: aa[16 * i + 4 * j + k] for i, x in enumerate(b) for j, y in enumerate(b) for k, z in enumerate(b)}
    seqs, name = {}, None
    for l in open(fa):
        l = l.strip()
        if l.startswith(">"): name = l[1:].split()[0]; seqs[name] = []
        elif name: seqs[name].append(l)
    rows = []
    with open(faa, "w") as out:
        for g, s in seqs.items():
            s = "".join(s).upper(); flags = []
            if len(s) % 3: flags.append("len_not_mult3")
            p = "".join(code.get(s[k:k + 3], "X") for k in range(0, len(s) - len(s) % 3, 3))
            if s[:3] in {"ATG", "GTG", "TTG"}: p = "M" + p[1:]
            else: flags.append("no_start")
            if p.endswith("*"): p = p[:-1]
            else: flags.append("no_stop")
            if "*" in p: flags.append(f"internal_stop:{p.count('*')}")
            if len(p) < 100: flags.append("short_lt100aa")
            rows.append({"Gene": g, "nt_len": len(s), "aa_len": len(p), "translation_flags": ";".join(flags)})
            out.write(f">{g}\n{p.replace('*', 'X')}\n")
    pd.DataFrame(rows).to_csv(qc, sep="\t", index=False)

# ---------------------------------------------------------------- parsers
def uniref_hits(path="res/uniref90.tsv"):
    cols = "qseqid sseqid pident length qlen slen evalue bitscore qcovhsp scovhsp stitle".split()
    ur = pd.read_csv(path, sep="\t", header=None, names=cols).sort_values(["qseqid", "bitscore"], ascending=[True, False])
    ur["rep"] = ur.sseqid.str.replace("UniRef90_", "", regex=False)
    best = ur.groupby("qseqid").head(1).copy()
    best["uniref_name"] = best.stitle.str.extract(r"^\S+ (.+?) n=\d")[0]
    best["uniref_tax"] = best.stitle.str.extract(r"Tax=(.+?) TaxID")[0]
    best = best.rename(columns={"qseqid": "Gene", "sseqid": "uniref90_best", "pident": "uniref90_pid",
                                "qcovhsp": "uniref90_qcov", "evalue": "uniref90_evalue"})
    # UniParc (UPI...) representatives carry no annotation: take xrefs from the best
    # UniProtKB-represented cluster among the top 5 hits instead
    xs = (ur[~ur.rep.str.startswith("UPI")].groupby("qseqid").head(1)[["qseqid", "rep", "pident"]]
          .rename(columns={"qseqid": "Gene", "rep": "up_acc", "pident": "xref_hit_pid"}))
    keep = ["Gene", "uniref90_best", "uniref90_pid", "uniref90_qcov", "uniref90_evalue", "uniref_name", "uniref_tax"]
    return best[keep].merge(xs, on="Gene", how="left")

def uniprot_xrefs(accs, chunk=200):
    fields = "accession,protein_name,go,xref_kegg,xref_eggnog,xref_pfam,xref_interpro,ec"
    s, out = requests.Session(), []
    for i in range(0, len(accs), chunk):
        for k in range(5):
            try:
                r = s.get("https://rest.uniprot.org/uniprotkb/accessions",
                          params={"accessions": ",".join(accs[i:i + chunk]), "fields": fields, "format": "tsv"}, timeout=60)
                if r.ok: out.append(pd.read_csv(io.StringIO(r.text), sep="\t")); break
                print("UniProt HTTP", r.status_code, r.text[:200], file=sys.stderr)
            except requests.RequestException as e:
                print("UniProt retry:", e, file=sys.stderr)
            time.sleep(3 * (k + 1))
    up = pd.concat(out, ignore_index=True)
    return up.rename(columns={"Entry": "up_acc", "Protein names": "up_protein", "Gene Ontology (GO)": "GO_uniprot",
                              "KEGG": "up_kegg_gene", "eggNOG": "up_eggnog", "Pfam": "up_pfam",
                              "InterPro": "up_interpro", "EC number": "up_EC"})

def cog_hits(path="res/cog.tsv"):
    cg = pd.read_csv(path, sep="\t", header=None,
                     names="qseqid sseqid pident length qlen slen evalue bitscore qcovhsp sstart send".split())
    cg["prot"] = cg.sseqid.str.replace(r"_(\d+)$", r".\1", regex=True)   # FASTA uses _1 for .1
    cc = pd.read_csv("db/cog-20.cog.csv", header=None, usecols=[2, 4, 6], names=["prot", "fp", "COG"], encoding="latin-1")
    m = cg.merge(cc[cc.prot.isin(set(cg.prot))], on="prot").dropna(subset=["COG"])
    def overlap(r):
        best = 0
        for seg in str(r.fp).split("="):
            try: a, b = map(int, seg.split("-")); best = max(best, min(b, r.send) - max(a, r.sstart))
            except ValueError: pass
        return best
    m["keep"] = m.apply(overlap, axis=1) >= 0.3 * m.length
    m = m[m.keep | ~m.qseqid.isin(m.loc[m.keep, "qseqid"])]
    cdef = pd.read_csv("db/cog-20.def.tab", sep="\t", header=None, encoding="latin-1", usecols=[0, 1, 2],
                       names=["COG", "COG_cat", "COG_name"])
    fd = pd.read_csv("db/fun-20.tab", sep="\t", header=None, names=["c", "col", "d"]).set_index("c").d.to_dict()
    j = lambda s, sep=";": sep.join(dict.fromkeys(s))
    g = (m.merge(cdef, on="COG", how="left").groupby("qseqid")
         .agg(COG=("COG", j), COG_cat=("COG_cat", j), COG_name=("COG_name", lambda s: j(s, "; ")),
              COG_pid=("pident", "first")).reset_index().rename(columns={"qseqid": "Gene"}))
    g["COG_cat_desc"] = g.COG_cat.apply(lambda s: "; ".join(dict.fromkeys(fd.get(c, "?") for c in "".join(s.split(";")))))
    return g

def pfam_hits(path="res/pfam.domtbl"):
    rows = []
    for l in open(path):
        if l.startswith("#"): continue
        f = l.split()
        rows.append({"Gene": f[3], "acc": f[1].split(".")[0], "name": f[0], "from": int(f[19]), "desc": " ".join(f[22:])})
    pf = pd.DataFrame(rows).sort_values(["Gene", "from"])
    p2g = {}
    for l in open("db/pfam2go"):
        m = re.match(r"Pfam:(PF\d+) \S+ > GO:(.+) ; (GO:\d+)", l.strip())
        if m: p2g.setdefault(m.group(1), []).append(f"{m.group(3)} {m.group(2)}")
    j = lambda s, sep=";": sep.join(dict.fromkeys(s))
    g = pf.groupby("Gene").agg(Pfam=("acc", j), Pfam_names=("name", j), Pfam_desc=("desc", lambda s: j(s, "; "))).reset_index()
    g["GO_pfam2go"] = g.Pfam.apply(lambda s: "; ".join(dict.fromkeys(x for a in s.split(";") for x in p2g.get(a, []))) or None)
    return g

def kofam_hits(path="res/kofam_detail.tsv", kegg=None):
    rows = []
    for l in open(path):
        if not l.startswith("*"): continue            # '*' = above KO-specific threshold
        f = l.rstrip("\n").split("\t")
        rows.append({"Gene": f[1].strip(), "KO": f[2].strip(), "ko_score": float(f[4]), "ko_evalue": float(f[5]),
                     "KO_definition": f[6].strip().strip("\"")})
    k = (pd.DataFrame(rows).sort_values(["Gene", "ko_score"], ascending=[True, False])
         .groupby("Gene").head(1).reset_index(drop=True))
    if kegg is not None:
        kg = pd.read_csv(kegg, sep="\t", dtype=str)
        pw = (kg.dropna(subset=["KO"]).groupby("KO")
              .agg(KEGG_level2=("Level2", lambda s: "; ".join(dict.fromkeys(s))),
                   KEGG_pathway=("Level3", lambda s: "; ".join(dict.fromkeys(s)))).reset_index())
        k = k.merge(pw, on="KO", how="left")
    return k

# ---------------------------------------------------------------- merge
def merge(a):
    qc = pd.read_csv("translation_qc.tsv", sep="\t")
    rt = pd.read_csv(a.rtab, sep="\t").set_index("Gene"); rt = rt[rt.sum(axis=1) > 0]
    meta = pd.read_csv(a.meta, sep="\t")
    cls = lambda n, N: np.select([n == N, n >= 2, n == 1], ["Core", "Accessory", "Unique"], "Absent")
    pc = pd.DataFrame({"Gene": rt.index, "n_isolates": rt.sum(axis=1).values}); pc["class_all"] = cls(pc.n_isolates, rt.shape[1])
    for pond in sorted(meta.Pond.unique()):
        lab = pond.replace(" ", ""); iso = meta.Group[meta.Pond == pond]
        n = rt[iso].sum(axis=1).values; pc[f"n_{lab}"] = n; pc[f"class_{lab}"] = cls(n, len(iso))
    sp = pd.read_csv(a.swissprot, sep="\t")[["qseqid", "stitle", "pident"]]
    sp = sp.assign(swissprot_best=sp.stitle.str.extract(r"Full=([^;\[]+)")[0].str.strip()) \
           .rename(columns={"qseqid": "Gene", "pident": "swissprot_pid"})[["Gene", "swissprot_best", "swissprot_pid"]]
    ur = uniref_hits()
    up = uniprot_xrefs(ur.up_acc.dropna().unique().tolist())
    A = (qc.merge(pc, on="Gene", how="left").merge(sp, on="Gene", how="left").merge(ur, on="Gene", how="left")
           .merge(up, on="up_acc", how="left").merge(kofam_hits(kegg=a.kegg), on="Gene", how="left")
           .merge(cog_hits(), on="Gene", how="left").merge(pfam_hits(), on="Gene", how="left"))
    # one best-available name per gene, with its source
    name_src = [("swissprot_best", "SwissProt"), ("KO_definition", "KOfam"), ("COG_name", "COG"),
                ("up_protein", "UniProtKB"), ("Pfam_desc", "Pfam")]
    bad = re.compile(r"^(uncharacteri[sz]ed|hypothetical|putative uncharacterized|DUF\d+|Domain of unknown function)", re.I)
    def pick(r):
        for c, s in name_src:
            v = r[c]
            if isinstance(v, str) and v and not bad.match(v): return pd.Series([v, s])
        return pd.Series([None, None])
    A[["best_name", "best_name_source"]] = A.apply(pick, axis=1)
    A.to_csv(a.out, sep="\t", index=False)
    cov = {c: int(A[c].notna().sum()) for c in ["swissprot_best", "uniref90_best", "GO_uniprot", "KO", "COG", "Pfam", "GO_pfam2go", "best_name"]}
    print(f"{len(A)} genes |", cov)
    return A

if __name__ == "__main__":
    if sys.argv[1] == "translate": translate(*sys.argv[2:5]); sys.exit()
    p = argparse.ArgumentParser(); p.add_argument("cmd")
    for o in ["--rtab", "--meta", "--kegg", "--swissprot", "--out"]: p.add_argument(o, required=True)
    merge(p.parse_args())
