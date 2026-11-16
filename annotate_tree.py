"""Copy BEAST-style [&key=value,...] node annotations from a NEXUS tree string onto
an already-parsed (ete3-style) tree, as node attributes (node.orientation, etc.).

brokilon's sranges code reads node.orientation, but read_nexus_trees doesn't keep the
metadata. Nodes are matched in post-order, which is the order they appear in the
newick string, so the parsed tree must come from the same string.
"""
import re
from pathlib import Path


def _split_meta(meta):
    """'a=1,b={x,y},c=z' -> {'a': '1', 'b': '{x,y}', 'c': 'z'} (commas inside {} kept)."""
    out, depth, buf = {}, 0, ""
    for ch in meta + ",":
        if ch == "{":
            depth += 1
        elif ch == "}":
            depth -= 1
        if ch == "," and depth == 0:
            if "=" in buf:
                k, v = buf.split("=", 1)
                out[k.strip()] = v.strip()
            buf = ""
        else:
            buf += ch
    return out


def postorder_annotations(newick):
    """List of annotation dicts, one per node, in post-order (newick text order of node ends)."""
    anns, i, n = [], 0, len(newick)
    while i < n:
        c = newick[i]
        if c in "(,":
            i += 1
            continue
        if c == ";":
            break
        # Node end: either a leaf label, or ')' followed by an optional internal label.
        is_leaf = c != ")"
        if c == ")":
            i += 1
        start = i
        while i < n and newick[i] not in "[:,();":
            i += 1
        label = newick[start:i]
        meta = {"__label__": label, "__leaf__": is_leaf}
        if i < n and newick[i] == "[":
            j = newick.index("]", i)
            body = newick[i + 1:j]
            meta.update(_split_meta(body[1:] if body.startswith("&") else body))
            i = j + 1
        if i < n and newick[i] == ":":
            i += 1
            while i < n and newick[i] not in ",();[":
                i += 1
            if i < n and newick[i] == "[":  # metadata after the branch length
                j = newick.index("]", i)
                body = newick[i + 1:j]
                meta.update(_split_meta(body[1:] if body.startswith("&") else body))
                i = j + 1
        anns.append(meta)
    return anns


def read_tree_string(path, index=0):
    text = Path(path).read_text()
    trees = re.findall(r"^\s*tree\s+\S+\s*=\s*(?:\[&[RU]\]\s*)?(.*;)", text, re.M | re.I)
    return trees[index]


def _species(name):
    return re.sub(r"_(first|last)$", "", name) if name else name


def _set(node, attr, v):
    if hasattr(node, "add_feature"):
        node.add_feature(attr, v)   # ete3
    else:
        setattr(node, attr, v)


def annotate(tree, path, index=0, keys=None, taxon_map=None):
    """Set node attributes from the file's annotations. Returns the tree.

    The CCD-output format (ancestral_range=X_first) is mapped to the BEAST format
    (range=X, taxon=X on leaves) so both trees look the same to brokilon.
    """
    anns = postorder_annotations(read_tree_string(path, index))
    nodes = list(tree.traverse("postorder"))
    if len(nodes) != len(anns):
        raise ValueError(f"{path}: parsed tree has {len(nodes)} nodes, string has {len(anns)}")
    for node, meta in zip(nodes, anns):
        if node.is_leaf() != meta.pop("__leaf__") or (node.is_leaf() and str(node.name) != meta["__label__"]):
            raise ValueError(f"{path}: node order mismatch at {node.name!r} vs {meta['__label__']!r}")
        meta.pop("__label__")
        for k, v in meta.items():
            if keys is None or k in keys:
                attr = re.sub(r"\W", "_", k)  # e.g. height_95%_HPD -> height_95__HPD
                _set(node, attr, v)
        # Normalise to the BEAST/SRFBD attribute names. brokilon tests these with
        # hasattr(), so only set them where they really exist -- never as None.
        rng = meta.get("range") or _species(meta.get("ancestral_range"))
        if rng:
            _set(node, "range", rng)
        taxon = meta.get("taxon")
        if taxon is None and node.is_leaf() and taxon_map is not None:
            taxon = _species(taxon_map[int(node.name)])
        if taxon:
            _set(node, "taxon", taxon)
    return tree