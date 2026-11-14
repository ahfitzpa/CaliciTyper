"""Configuration for the recombination module.

YAML or TOML on disk, a dataclass in memory, with every threshold, weight and
tool path overridable. Defaults are the ones documented in the README; nothing
here is a magic number buried in the code.
"""
from __future__ import annotations

import os
import shutil
from dataclasses import dataclass, field, asdict
from pathlib import Path

ALL_METHODS = ("rdp", "geneconv", "bootscan", "maxchi", "chimaera", "siscan", "threeseq")

# Breakpoint regions, per virus. These say only whether a breakpoint falls
# somewhere biologically interesting; boundaries are strain-dependent and the
# windows are deliberately generous.
#
# Norovirus: ORF1 (nonstructural, including RdRp), ORF2 (VP1 capsid), ORF3
# (VP2). The dominant recombination hotspot is the ORF1/ORF2 junction.
#
# Sapovirus is NOT the same shape and cannot reuse the norovirus model. ORF1
# encodes the nonstructural proteins AND the major capsid VP1 as one
# polyprotein; ORF2 encodes the minor capsid VP2. The hotspot is therefore the
# NS7(RdRp)/VP1 junction, which lies INSIDE ORF1 rather than between two ORFs
# -- a porcine GIII recombinant maps its breakpoint to nt 5139 at exactly that
# boundary. Labelling that an "ORF1-ORF2 junction" would be wrong, even though
# the coordinate happens to be similar.
REGIONS_NOROVIRUS = {
    "ORF1-ORF2 junction": [4800, 5600],
    "ORF1 (nonstructural)": [5, 5100],
    "ORF2 (VP1 capsid)": [5000, 6800],
    "ORF3 (VP2)": [6700, 7600],
}

REGIONS_SAPOVIRUS = {
    "RdRp-VP1 junction": [4800, 5500],
    "ORF1 nonstructural (NS1-NS7)": [10, 5100],
    "ORF1 VP1 capsid": [5100, 6800],
    "ORF2 (VP2)": [6700, 7200],
}

REGIONS_BY_VIRUS = {"norovirus": REGIONS_NOROVIRUS, "sapovirus": REGIONS_SAPOVIRUS}
DEFAULT_REGIONS = REGIONS_NOROVIRUS


@dataclass
class Scoring:
    method_support: float = 0.35
    phylogenetic_support: float = 0.20
    breakpoint_consistency: float = 0.20
    parental_distinctness: float = 0.15
    biological_plausibility: float = 0.10

    def weights(self):
        return asdict(self)

    def validate(self):
        total = sum(self.weights().values())
        if abs(total - 1.0) > 1e-6:
            raise ValueError(f"scoring weights must sum to 1.0, got {total:.4f}")


@dataclass
class StatusRules:
    weak_min_methods: int = 1
    possible_min_methods: int = 2
    likely_min_methods: int = 3
    likely_require_phylo: bool = True
    strong_min_methods: int = 3
    strong_require_phylo: bool = True
    strong_require_distinct: bool = True


@dataclass
class RecombConfig:
    min_sequence_length: int = 5000
    min_informative_sites: int = 500
    max_ambiguous_fraction: float = 0.05
    pvalue_threshold: float = 0.05
    multiple_testing_correction: str = "bonferroni"       # bonferroni | none
    max_breakpoint_distance: int = 250
    methods: tuple = ALL_METHODS
    run_gard: bool = False                                 # stage 3; off until validated
    build_segment_trees: bool = True
    segment_bootstraps: int = 100
    min_segment_sites: int = 100
    # Columns either side of a breakpoint used for the phylogeny check.
    # 0 means "use the whole segment", which confirms recombination but
    # cannot localise it -- see segment_support().
    segment_window: int = 1000
    seed: int = 3
    virus: str = "norovirus"              # norovirus | sapovirus
    regions: dict = field(default_factory=dict)   # empty -> taken from `virus`
    scoring: Scoring = field(default_factory=Scoring)
    status_rules: StatusRules = field(default_factory=StatusRules)
    openrdp_exe: str = ""                                  # else $CALICITYPER_OPENRDP, else PATH
    hyphy_exe: str = ""

    def validate(self):
        if self.virus not in REGIONS_BY_VIRUS:
            raise ValueError(f"virus must be one of {sorted(REGIONS_BY_VIRUS)}, got {self.virus!r}")
        if not self.regions:
            self.regions = dict(REGIONS_BY_VIRUS[self.virus])
        self.scoring.validate()
        bad = [m for m in self.methods if m not in ALL_METHODS]
        if bad:
            raise ValueError(f"unknown method(s): {bad}; permitted: {list(ALL_METHODS)}")
        if self.multiple_testing_correction not in ("bonferroni", "none"):
            raise ValueError("multiple_testing_correction must be 'bonferroni' or 'none'")
        if not 0 < self.pvalue_threshold <= 1:
            raise ValueError("pvalue_threshold must be in (0, 1]")
        return self

    def openrdp(self):
        """Path to the OpenRDP CLI, or None.

        OpenRDP pins numpy<2 while calicityper needs numpy>=2, so it cannot share
        this environment. It is therefore treated like BLAST and MAFFT: an
        external executable, typically in its own conda env, located by
        $CALICITYPER_OPENRDP or PATH.
        """
        cand = self.openrdp_exe or os.environ.get("CALICITYPER_OPENRDP") or shutil.which("openrdp")
        return cand if cand and Path(cand).exists() else None

    def hyphy(self):
        cand = self.hyphy_exe or os.environ.get("CALICITYPER_HYPHY") or shutil.which("hyphy")
        return cand if cand and Path(cand).exists() else None


def _flatten_load(data):
    """Map the nested on-disk layout onto RecombConfig."""
    cfg = RecombConfig()
    rec = (data or {}).get("recombination", {}) or {}
    for k, v in rec.items():
        if k == "methods":
            cfg.methods = tuple(str(x).lower() for x in v)
        elif hasattr(cfg, k):
            setattr(cfg, k, v)
    sc = (data or {}).get("scoring", {}) or {}
    for k, v in sc.items():
        if hasattr(cfg.scoring, k):
            setattr(cfg.scoring, k, float(v))
    sr = (data or {}).get("status_rules", {}) or {}
    for status, rules in sr.items():
        for k, v in (rules or {}).items():
            attr = {"min_methods": f"{status}_min_methods",
                    "require_phylogenetic_support": f"{status}_require_phylo",
                    "require_parental_distinctness": f"{status}_require_distinct"}.get(k)
            if attr and hasattr(cfg.status_rules, attr):
                setattr(cfg.status_rules, attr, v)
    if "regions" in rec:
        cfg.regions = {str(k): list(v) for k, v in rec["regions"].items()}
    return cfg.validate()


def load(path=None):
    """Load config from YAML or TOML; defaults when path is None."""
    if not path:
        return RecombConfig().validate()
    p = Path(path)
    if not p.exists():
        raise FileNotFoundError(f"config not found: {p}")
    text = p.read_text()
    if p.suffix.lower() in (".toml",):
        import tomllib
        data = tomllib.loads(text)
    else:
        import yaml
        data = yaml.safe_load(text)
    return _flatten_load(data)
