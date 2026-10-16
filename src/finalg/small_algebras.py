# ===============================================
#   Small algebras: a catalogue of small groups
# ===============================================
"""A catalogue of all groups of order 1 through 40, one for each isomorphism class (181 groups in all).

The groups are stored explicitly, as pretty-printed JSON files (one per group, plus a catalogue,
``index.json``), in the package ``finalg.data.small_groups``. They are located with
``importlib.resources``, so this works for a source checkout, an editable install, and an installed
wheel alike. How each group can be built from simpler ones is shown, with finalg code, in the notebooks
in ``notebooks/small_groups`` (three groups of order 32 cannot be built that way, and were converted
from GAP; see ``tools/gap_to_finalg.g``).

The functions are:

* ``small_groups(order=None)`` returns the dictionary of all small groups, keyed by order, or, if an order
  is given, the list of all groups of that order. The groups are created once and shared.
* ``identify_small_group(alg)`` returns the small group that a given group is isomorphic to, or None.
* ``load_group(key)`` loads one group from its JSON file, by name, by SmallGroup id, or by file path, and
  ``load_groups(...)`` loads many of them (all of them by default), as a list or a dictionary.
* ``groups_index()`` returns the catalogue, ``index.json``, as a list of dictionaries, ``data_dir()``
  returns the folder that holds the data, and ``verify_library()`` checks all of it.

Group identifiers follow GAP's SmallGroup library: ``(order, index)``, e.g. ``(8, 4)`` is the
quaternion group Q8, and within an order the groups are listed in the order of GAP's library. File names
look like ``G08_04_Q8.json``.

Group names follow the conventions of the generators: ``Z6`` is cyclic of order 6, ``D4`` is the
dihedral group of order 8 (the symmetries of a square), ``Q8``, ``Q16`` and ``Q32`` are the quaternion
and generalized quaternion groups, ``Dic3`` and ``Dic5`` are the dicyclic groups of orders 12 and 20,
``Z4_x_Z2`` is a direct product, ``Z7_sdp_Z3`` is a semi-direct product, ``S4`` is the symmetric group on
4 letters, and so on. The three groups of order 32 that were converted from GAP are named ``SG32_8``,
``SG32_15`` and ``SG32_32``.
"""

import json
import re
from functools import lru_cache
from importlib import resources
from pathlib import Path

from finalg.group import Group
from finalg.make_finite_algebra import make_finite_algebra

# The package that holds index.json and the json folder
DATA_PACKAGE = "finalg.data.small_groups"

MAX_SMALL_ORDER = 40

# The known number of groups of order 1, 2, ..., 40 (OEIS A000001), used as a sanity check.
_GROUP_COUNTS = (1, 1, 1, 2, 1, 2, 1, 5, 2, 2, 1, 5, 1, 2, 1, 14, 1, 5, 1, 5,
                 2, 2, 1, 15, 2, 2, 5, 4, 1, 4, 1, 51, 1, 2, 1, 14, 1, 2, 2, 14)


# ----------------------------------------------------------------------------------------
#  Reading the data
# ----------------------------------------------------------------------------------------

_index_cache = {}


def data_dir(directory=None):
    """Return the folder that holds index.json and the json folder of the library: the package
    ``finalg.data.small_groups`` or, if `directory` is given, that directory instead. The result
    supports ``/`` and ``read_text()`` (it is a ``Path`` for a source checkout or a regular install)."""
    return Path(directory) if directory else resources.files(DATA_PACKAGE)


def _read_json(directory, *parts):
    """Parse the JSON file at data_dir(directory)/parts. Reading the text, instead of opening a file
    by path, keeps this independent of where, or in what form, the package is installed."""
    node = data_dir(directory)
    for part in parts:
        node = node / part
    return json.loads(node.read_text(encoding='utf-8'))


def _load_file(directory, filename):
    """Create the algebra that is stored in the library file `filename`."""
    return make_finite_algebra(_read_json(directory, 'json', filename))


def groups_index(directory=None):
    """Return the catalogue of the library as a list of dictionaries, one per group, in order of
    (order, SmallGroup index). Keys: order, gap_id, name, file, structure, abelian, source."""
    if directory:      # a user-supplied folder may change between calls, so it is never cached
        return [dict(r) for r in _read_json(directory, 'index.json')]
    if 'default' not in _index_cache:
        _index_cache['default'] = _read_json(None, 'index.json')
    return [dict(r) for r in _index_cache['default']]


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
    Group (or other algebra). Every call reads the file again and returns a new object; use
    ``small_groups`` for groups that are created once and shared."""
    if isinstance(key, Path) or (isinstance(key, str) and key.endswith('.json') and Path(key).exists()):
        return make_finite_algebra(str(key))
    record = _resolve(key, groups_index(directory))
    return _load_file(directory, record['file'])


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
    loaded = [((r['order'], r['gap_id']), _load_file(directory, r['file'])) for r in records]
    return dict(loaded) if as_dict else [g for _, g in loaded]


def verify_library(directory=None, verbose=True):
    """Load every group in the library and check it against index.json: the file exists, the algebra is a
    Group of the stated order and name, names and ids are unique, and there are exactly as many groups of
    each order as there should be (OEIS A000001). Raises AssertionError on a problem; returns the number of
    groups checked. (This does not test that the groups are pairwise non-isomorphic. The library was
    checked against GAP's SmallGroup ids when it was built.)"""
    index = groups_index(directory)
    assert len({r['name'] for r in index}) == len(index), "group names are not unique"
    assert len({(r['order'], r['gap_id']) for r in index}) == len(index), "SmallGroup ids are not unique"
    for r in index:
        alg = _load_file(directory, r['file'])
        assert type(alg).__name__ == 'Group', f"{r['file']} is a {type(alg).__name__}, not a Group"
        assert alg.order == r['order'], f"{r['file']} has order {alg.order}, expected {r['order']}"
        assert alg.name == r['name'], f"{r['file']} is named {alg.name!r}, expected {r['name']!r}"
    for n in sorted({r['order'] for r in index}):
        count = sum(1 for r in index if r['order'] == n)
        assert count == _GROUP_COUNTS[n - 1], f"order {n}: {count} groups, expected {_GROUP_COUNTS[n - 1]}"
    if verbose:
        print(f"OK: {len(index)} groups of orders {index[0]['order']} through {index[-1]['order']}.")
    return len(index)


# ----------------------------------------------------------------------------------------
#  The catalogue of small groups
# ----------------------------------------------------------------------------------------

@lru_cache(maxsize=None)
def _catalogue():
    """Load, once, the dictionary that maps each order from 1 to 40 to the list of all groups of
    that order, up to isomorphism, in the order of GAP's SmallGroup library."""
    groups = {n: [] for n in range(1, MAX_SMALL_ORDER + 1)}
    for r in sorted(groups_index(), key=lambda rec: (rec['order'], rec['gap_id'])):
        groups[r['order']].append(_load_file(None, r['file']))
    for n, grps in groups.items():
        if len(grps) != _GROUP_COUNTS[n - 1] or any(g.order != n for g in grps):
            raise RuntimeError(f"The catalogue of small groups of order {n} is wrong: "
                               f"{len(grps)} groups, but there should be {_GROUP_COUNTS[n - 1]}.")
    return groups


def small_groups(order=None):
    """Return the small groups, which are all the groups of order 1 through 40, up to isomorphism.

    If `order` is None (the default), a dictionary is returned whose keys are the orders 1, ..., 40 and
    whose values are the lists of all groups of that order, e.g., small_groups()[4] is [Z4, V4].
    If `order` is an integer from 1 to 40, then just the list of groups of that order is returned,
    e.g., small_groups(4) is [Z4, V4]. Within an order, the groups are in the order of GAP's SmallGroup
    library, so that the group at position i in the list has the SmallGroup id (order, i + 1).

    The Group objects are created once and then shared, so a group returned by one call is the very
    same object returned by another. The dictionary and lists are new each time, so changing them does
    not change the catalogue."""
    catalogue = _catalogue()
    if order is None:
        return {n: list(grps) for n, grps in catalogue.items()}
    if isinstance(order, bool) or not isinstance(order, int) or not 1 <= order <= MAX_SMALL_ORDER:
        raise ValueError(f"order = {order!r} is invalid; it must be an integer from 1 to {MAX_SMALL_ORDER}, "
                         f"or None for all of the small groups.")
    return list(catalogue[order])


def identify_small_group(alg):
    """Return the small group that the given algebra is isomorphic to, or None if there isn't one.

    There is one only if `alg` is a Group of order 40 or less, and then there is exactly one, because
    the small groups contain every isomorphism class once. The result is one of the shared objects
    returned by ``small_groups``. Isomorphism is tested with `fast_isomorphic`, which first compares
    invariants of the elements, so that, in practice, there are few, if any, groups of the same order
    that it has to try."""
    if not isinstance(alg, Group) or alg.order > MAX_SMALL_ORDER:
        return None
    candidates = _catalogue()[alg.order]
    for small in candidates:
        if alg is small:
            return small
    for small in candidates:
        if alg.fast_isomorphic(small):
            return small
    return None
