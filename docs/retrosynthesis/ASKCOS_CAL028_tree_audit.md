# ASKCOS CAL-028 retrosynthesis: source and audit

Source: [ASKCOS shared result](https://askcos.mit.edu/results?shared=0a31ab2d-cd9b-4b72-a129-7f43a8fd4bfb).

Unmodified 200-tree graph export: [ASKCOS_CAL028_treeResults.json](ASKCOS_CAL028_treeResults.json).

- SHA-256 of source JSON: `43b7d9f6740a407403c0087d361e3dd67f1e1bb2b5c18fd9a699c46a24aeed74`.
- Number of graph objects: **200**.
- `graph.num_reactions`: **56** with 4 steps; **144** with 5 steps.
- Number of distinct reaction-SMILES among all `nodes` of type `reaction`: **22**.
- Reaction adjacent to CAL-028 target in each tree: **68** proposals from methyl iodide (`CI`), **66** from formaldehyde (`C=O`), and **66** from methanol (`CO`).
- A tree is one predicted retrosynthetic route. Recurrent transformations across multiple trees are not counted as independent reaction types.

The independently proposed chemically convergent route is described in [CAL028_proposed_scheme.md](CAL028_proposed_scheme.md). Its choice of N-methylated stereodefined amino-acid building block avoids relying on the target-stage methylation that dominates ASKCOS recommendations.

## Count verification (Python)
```python
import json
from collections import Counter
from pathlib import Path
trees = json.loads(Path("ASKCOS_CAL028_treeResults.json").read_text())
assert len(trees) == 200
steps = Counter(t["graph"]["num_reactions"] for t in trees)
assert steps == {4: 56, 5: 144}
unique_rx = {node["smiles"] for t in trees for node in t["nodes"] if node["type"] == "reaction"}
assert len(unique_rx) == 22
root = "00000000-0000-0000-0000-000000000000"
reagents = Counter()
for t in trees:
    edge = next(e for e in t["edges"] if e["from"] == root)
    step = next(n["smiles"] for n in t["nodes"] if n["id"] == edge["to"])
    reagents[step.split(">>")[0].split(".")[0]] += 1
assert reagents == {"CI": 68, "C=O": 66, "CO": 66}
print("200 trees; 56 four-step; 144 five-step; 22 reactions; 68/66/66 reagents")
```
