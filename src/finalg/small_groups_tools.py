"""
small_groups_tools: helpers for the finalg small-groups library (all groups of order 1 through 40).

Three things live here:

1. Construction helpers used by the notebooks:
   * ``mult_action(modulus, k)``            -- action of Z_m on Z_modulus by multiplication
   * ``action_from_generators(N, H, images)`` -- an action of H on N by automorphisms, given the
                                                 images of generators (for ``Group.semidirect_product``)

2. Pretty-printed JSON output:
   * ``pretty_json(d)`` and ``dump_pretty(alg, path)`` write an algebra in the same row-per-line layout
     as the example files in ``finalg/data/algebras``.

3. Loading the serialized groups in the ``json`` folder next to this file:
   * ``load_group(key)``   -- one group, by name, by SmallGroup id (order, index), or by file path
   * ``load_groups(...)``  -- many groups (all of them by default) as a list or a dictionary
   * ``groups_index()``    -- the catalogue (index.json) as a list of dictionaries
   * ``verify_library()``  -- load everything and check it against the index

Group identifiers follow GAP's SmallGroup library: ``(order, index)``, e.g. ``(8, 4)`` is the
quaternion group Q8. File names look like ``G08_04_Q8.json``.
"""

import json
import re
from pathlib import Path

from finalg import make_finite_algebra

HERE = Path(__file__).resolve().parent
JSON_DIR = HERE / 'json'
INDEX_FILE = HERE / 'index.json'


# ----------------------------------------------------------------------------------------
#  Construction helpers
# ----------------------------------------------------------------------------------------

def mult_action(modulus, k):
    """The action of Z_m on Z_modulus (elements are strings of integers) in which h acts by
    x -> k^h * x (mod modulus). Use with ``Group.semidirect_product``."""
    return lambda h, x: str((int(x) * pow(k, int(h), modulus)) % modulus)


def _extend_hom(N, gen_images):
    """Extend {generator: image} to a map on all of N by f(a*g) = f(a)*f(g)."""
    if not gen_images:
        return {n: n for n in N.elements}                 # no images given: the identity automorphism
    f = {N.identity: N.identity}
    f.update(gen_images)
    queue = list(f)
    while queue:
        a = queue.pop()
        for g, img in gen_images.items():
            p = N.op(a, g)
            if p not in f:
                f[p] = N.op(f[a], img)
                queue.append(p)
    return f


def action_from_generators(N, H, images):
    """Build an action of the group H on the group N by automorphisms, from generators.

    `images` maps each of a set of generators of H to a dictionary {n_gen: image, ...} that says where the
    corresponding automorphism of N sends a set of generators of N (an empty dictionary means that the
    H-generator acts trivially). The result is a dictionary {h: {n: phi_h(n)}} that covers every element
    h of H and n of N, ready to be passed to ``N.semidirect_product(H, action)``, which verifies that it
    really is a homomorphism from H into the automorphisms of N."""
    identity_map = {n: n for n in N.elements}
    gen_autos = {h: _extend_hom(N, imgs) for h, imgs in images.items()}
    phi = {H.identity: identity_map}
    queue = [H.identity]
    while queue:
        h = queue.pop()
        for g, auto in gen_autos.items():
            hg = H.op(h, g)
            if hg not in phi:
                phi[hg] = {n: phi[h][auto[n]] for n in N.elements}     # phi_(hg) = phi_h o phi_g
                queue.append(hg)
    if len(phi) != H.order:
        raise ValueError(f"The generators given for {H.name} generate only {len(phi)} of its {H.order} elements.")
    return phi


# ----------------------------------------------------------------------------------------
#  Pretty JSON
# ----------------------------------------------------------------------------------------

def _wrapped_list(items, indent, width):
    """JSON-encode a list of scalars on as few lines as fit within `width` columns."""
    pieces = [json.dumps(x, ensure_ascii=False) for x in items]
    one_line = '[' + ', '.join(pieces) + ']'
    if indent + len(one_line) <= width:
        return one_line
    lines, current = [], ''
    for p in pieces:
        candidate = (current + ', ' + p) if current else p
        if current and indent + 2 + len(candidate) + 1 > width:
            lines.append(current + ',')
            current = p
        else:
            current = candidate
    lines.append(current)
    pad = ' ' * (indent + 2)
    return '[\n' + '\n'.join(pad + ln for ln in lines) + '\n' + ' ' * indent + ']'


def pretty_json(d, width=100):
    """Return the JSON text for an algebra dictionary (from ``alg.to_dict()``), laid out like the example
    files in finalg/data/algebras: one key per line, the elements wrapped to `width` columns, and one
    row of each Cayley table per line."""
    parts = []
    for key, value in d.items():
        k = json.dumps(key)
        if key in ('table', 'table2'):
            rows = ',\n'.join('    ' + json.dumps(list(r)) for r in value)
            parts.append(f'  {k}: [\n{rows}\n  ]')
        elif key == 'elements':
            parts.append(f'  {k}: ' + _wrapped_list(list(value), 2, width))
        else:
            parts.append(f'  {k}: ' + json.dumps(value, ensure_ascii=False))
    return '{\n' + ',\n'.join(parts) + '\n}\n'


def dump_pretty(alg, path, width=100):
    """Write an algebra to `path` as pretty-printed JSON that ``make_finite_algebra(path)`` can read back."""
    with open(path, 'w', encoding='utf-8') as fout:
        fout.write(pretty_json(alg.to_dict(), width))


# ----------------------------------------------------------------------------------------
#  Loading
# ----------------------------------------------------------------------------------------

_index_cache = {}


def groups_index(directory=None):
    """Return the catalogue of the library as a list of dictionaries, one per group, in order of
    (order, SmallGroup index). Keys: order, gap_id, name, file, structure, abelian, source."""
    base = Path(directory) if directory else HERE
    key = str(base)
    if key not in _index_cache:
        with open(base / 'index.json', encoding='utf-8') as fin:
            _index_cache[key] = json.load(fin)
    return [dict(r) for r in _index_cache[key]]


def _json_dir(directory):
    return (Path(directory) / 'json') if directory else JSON_DIR


def _resolve(key, index):
    """Find the index record for a name, a (order, index) pair, or a 'G08_04' / 'G08_04_Q8.json' style key."""
    if isinstance(key, (tuple, list)) and len(key) == 2:
        order, gid = int(key[0]), int(key[1])
        for r in index:
            if r['order'] == order and r['gap_id'] == gid:
                return r
        raise KeyError(f"There is no group with SmallGroup id ({order}, {gid}) in the library.")
    if isinstance(key, str):
        for r in index:
            if r['name'] == key or r['file'] == key or r['file'][:-5] == key:
                return r
        m = re.fullmatch(r'G(\d+)_(\d+)', key)
        if m:
            return _resolve((int(m.group(1)), int(m.group(2))), index)
        close = [r['name'] for r in index if key.lower() in r['name'].lower()][:8]
        hint = f" Did you mean one of {close}?" if close else ''
        raise KeyError(f"There is no group named {key!r} in the library.{hint}")
    raise TypeError(f"A key must be a group name (str) or a (order, index) pair, not {key!r}.")


def load_group(key, directory=None):
    """Load a single group from the library.

    `key` can be a name ('Q8', 'F21', ...), a SmallGroup id as a pair (order, index) such as (8, 4),
    a file stem such as 'G08_04', or the path of any JSON file that finalg can read. Returns a finalg
    Group (or other algebra)."""
    if isinstance(key, Path) or (isinstance(key, str) and key.endswith('.json') and Path(key).exists()):
        return make_finite_algebra(str(key))
    record = _resolve(key, groups_index(directory))
    return make_finite_algebra(str(_json_dir(directory) / record['file']))


def load_groups(orders=None, keys=None, names=None, as_dict=False, directory=None):
    """Load many groups from the library at once.

    * With no arguments, all of the groups are loaded (orders 1 through 40; 181 groups).
    * `orders`: an order (int) or an iterable of orders, e.g. 16 or range(1, 21).
    * `names`: an iterable of group names, e.g. ['Q8', 'D4', 'Pauli'].
    * `keys`: an iterable of keys of any kind accepted by ``load_group``, e.g. [(8, 4), 'F21'].

    If both `orders` and `names` are given, only the groups satisfying both are loaded. `keys`, if given,
    selects exactly those groups, in that order.

    The result is a list of groups, in order of (order, SmallGroup index), unless `as_dict` is True, in which
    case it is a dictionary whose keys are the SmallGroup ids (order, index) and whose values are the groups."""
    index = groups_index(directory)
    if keys is not None:
        records = [_resolve(k, index) for k in keys]
    else:
        records = index
        if orders is not None:
            wanted = {orders} if isinstance(orders, int) else set(orders)
            records = [r for r in records if r['order'] in wanted]
        if names is not None:
            wanted_names = {names} if isinstance(names, str) else set(names)
            missing = wanted_names - {r['name'] for r in index}
            if missing:
                raise KeyError(f"There are no groups named {sorted(missing)} in the library.")
            records = [r for r in records if r['name'] in wanted_names]
    folder = _json_dir(directory)
    loaded = [((r['order'], r['gap_id']), make_finite_algebra(str(folder / r['file']))) for r in records]
    return dict(loaded) if as_dict else [g for _, g in loaded]


def verify_library(directory=None, verbose=True):
    """Load every group in the library and check it against index.json: the file exists, the algebra is a
    Group of the stated order and name, names and ids are unique, and there are exactly as many groups of
    each order as there should be (OEIS A000001). Raises AssertionError on a problem; returns the number of
    groups checked."""
    expected = [1, 1, 1, 2, 1, 2, 1, 5, 2, 2, 1, 5, 1, 2, 1, 14, 1, 5, 1, 5,
                2, 2, 1, 15, 2, 2, 5, 4, 1, 4, 1, 51, 1, 2, 1, 14, 1, 2, 2, 14]
    index = groups_index(directory)
    folder = _json_dir(directory)
    assert len({r['name'] for r in index}) == len(index), "group names are not unique"
    assert len({(r['order'], r['gap_id']) for r in index}) == len(index), "SmallGroup ids are not unique"
    for r in index:
        alg = make_finite_algebra(str(folder / r['file']))
        assert type(alg).__name__ == 'Group', f"{r['file']} is a {type(alg).__name__}, not a Group"
        assert alg.order == r['order'], f"{r['file']} has order {alg.order}, expected {r['order']}"
        assert alg.name == r['name'], f"{r['file']} is named {alg.name!r}, expected {r['name']!r}"
    for n in sorted({r['order'] for r in index}):
        count = sum(1 for r in index if r['order'] == n)
        assert count == expected[n - 1], f"order {n}: {count} groups, expected {expected[n - 1]}"
    if verbose:
        print(f"OK: {len(index)} groups of orders {index[0]['order']} through {index[-1]['order']}.")
    return len(index)
