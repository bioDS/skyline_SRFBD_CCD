"""Rewrite NEXUS tree files so leaves are labelled by integer IDs from ONE shared
Translate map, and return that map as {int: taxon_name} for brokilon's sranges code.

brokilon's prelabel_tree does `taxon_map[int(leaf.name)]`, so:
  * leaf names in the newick string must be integers, and
  * taxon_map must be a dict {int: name}, not a tree.
Both trees must use the SAME numbering, otherwise clades won't match.
"""
import re
from pathlib import Path

# A leaf label sits right after "(" or "," and right before "[" (metadata) or ":" (branch length).
LEAF_RE = re.compile(r"(?<=[(,])([A-Za-z][A-Za-z0-9_.\-]*)(?=[\[:,)])")


def read_taxlabels(path):
    text = Path(path).read_text()
    block = re.search(r"Taxlabels(.*?);", text, re.S | re.I).group(1)
    return block.split()


def build_taxon_map(path):
    """{1: 'Acynodon_iberoccitanus_first', 2: ...} from a file's Taxlabels."""
    return {i + 1: name for i, name in enumerate(read_taxlabels(path))}

SUFFIX_SWAP = {"_first": "_last", "_last": "_first"}


def renumber(in_path, out_path, taxon_map, fix_suffixes=True):
    name_to_id = {v: k for k, v in taxon_map.items()}
    text = Path(in_path).read_text()
    lines = text.splitlines()
    tree_re = re.compile(r"(\s*tree\s+\S+\s*=\s*)(.*)", re.I)

    # All leaf names actually used in the tree(s), so a swap can't merge two real ranges.
    all_leaves = {n for line in lines if (m := tree_re.match(line))
                  for n in LEAF_RE.findall(m.group(2))}

    def resolve(name):
        if name in name_to_id:
            return name
        if fix_suffixes:
            for suf, other in SUFFIX_SWAP.items():
                if name.endswith(suf):
                    alt = name[: -len(suf)] + other
                    if alt in name_to_id and alt not in all_leaves:
                        print(f"  {in_path}: renaming {name} -> {alt}")
                        return alt
        raise KeyError(f"{in_path}: leaf {name!r} not in taxon map")

    seen = set()

    def sub_leaves(newick):
        def repl(m):
            name = resolve(m.group(1))
            seen.add(name)
            return str(name_to_id[name])
        return LEAF_RE.sub(repl, newick)

    out_lines = []
    for line in lines:
        m = tree_re.match(line)
        out_lines.append(m.group(1) + sub_leaves(m.group(2)) if m else line)

    missing = set(name_to_id) - seen
    if missing:
        raise ValueError(f"{in_path}: {len(missing)} taxa never appear as leaves, e.g. {sorted(missing)[:3]}")

    translate = "\tTranslate\n" + ",\n".join(f"\t\t{k} {v}" for k, v in sorted(taxon_map.items())) + "\n;"
    out = re.sub(r"\tTranslate.*?\n;", lambda _: translate, "\n".join(out_lines) + "\n", count=1, flags=re.S)
    Path(out_path).write_text(out)
    return out_path

if __name__ == "__main__":
    import sys
    ref, *pairs = sys.argv[1:]
    tmap = build_taxon_map(ref)
    for i in range(0, len(pairs), 2):
        renumber(pairs[i], pairs[i + 1], tmap)
        print("wrote", pairs[i + 1])