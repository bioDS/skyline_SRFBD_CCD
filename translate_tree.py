"""Replace numeric taxon IDs in a Nexus tree file with names from its Translate block.

Usage:
    python translate_tree.py input.tre output.tre            # Nexus out, metadata kept
    python translate_tree.py input.tre output.nwk --newick   # plain Newick, metadata stripped
"""
import re
import sys
from pathlib import Path

# A tip ID: digits right after "(" or ",", followed by "[", ":", ",", ")" or end of chunk
TIP = re.compile(r"(?<=[(,])(\d+)(?=[\[:,)]|$)")
COMMENT = re.compile(r"(\[.*?\])", re.S)


def parse_translate(text):
    block = re.search(r"\bTranslate\s+(.*?);", text, re.S | re.I)
    if not block:
        raise ValueError("No Translate block found")
    mapping = {}
    for entry in block.group(1).split(","):
        parts = entry.split(None, 1)
        if len(parts) == 2:
            mapping[parts[0]] = parts[1].strip().strip("'\"")
    return mapping


def translate_newick(newick, mapping):
    # Split out [&...] comments so numbers inside metadata are never touched
    chunks = COMMENT.split(newick)
    for i in range(0, len(chunks), 2):
        chunks[i] = TIP.sub(lambda m: mapping.get(m.group(1), m.group(1)), chunks[i])
    return "".join(chunks)


def strip_metadata(newick):
    newick = COMMENT.sub("", newick)
    return re.sub(r"\)\d+", ")", newick)  # drop numeric internal node labels


def main(in_path, out_path, newick_only=False):
    text = Path(in_path).read_text()
    mapping = parse_translate(text)

    trees = re.findall(r"^\s*tree\s+(\S+)\s*=\s*(.*?;)\s*$", text, re.M | re.I)
    if not trees:
        raise ValueError("No tree statements found")

    if newick_only:
        out = "\n".join(strip_metadata(translate_newick(t, mapping)) for _, t in trees)
        Path(out_path).write_text(out + "\n")
    else:
        def replace_tree(m):
            return m.group(1) + translate_newick(m.group(2), mapping)

        text = re.sub(r"(^\s*tree\s+\S+\s*=\s*)(.*?;)", replace_tree, text, flags=re.M | re.I)
        text = re.sub(r"\s*\bTranslate\s+.*?;", "", text, count=1, flags=re.S | re.I)
        Path(out_path).write_text(text)

    print(f"Translated {len(trees)} tree(s) using {len(mapping)} taxa -> {out_path}")


if __name__ == "__main__":
    if len(sys.argv) < 3:
        sys.exit(__doc__)
    main(sys.argv[1], sys.argv[2], newick_only="--newick" in sys.argv[3:])