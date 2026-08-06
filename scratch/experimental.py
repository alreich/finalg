# ---------------------------------
#         EXPERIMENTAL
# Find all groups of a given order
# ---------------------------------

def _no_conflict(p1, p2):
    """Returns True only if no element of p1 equals the corresponding element of p2."""
    return all([p1[i] != p2[i] for i in range(len(p1))])


def _no_conflicts(items):
    """Return True if each possible pair, from a list of items, has no conflicts."""
    return all(_no_conflict(combo[0], combo[1]) for combo in it.combinations(items, 2))


def _filter_out_conflicts(perms, perm, n):
    """Filter out all permutations in perms that conflict with perm,
    and don't have n as the first element."""
    nperms = [q for q in perms if q[0] == n]
    return [p for p in nperms if _no_conflict(p, perm)]


def generate_all_group_tables(order):
    """Experimental Code: Return a list of all arrays that correspond to multiplication tables for groups
    of a specific order.

    WARNING: The algorithm here is inefficient, so even very small values of 'order' will result in very
    long runtimes (e.g., order=6 ==> ~5-6 hrs runtime).
    """
    row0 = list(range(order))
    row_candidates = [[row0]]
    for row_num in range(1, order):
        row_candidates.append(_filter_out_conflicts(it.permutations(row0), row0, row_num))
    table_candidates = [cand for cand in it.product(*row_candidates) if is_table_associative(list(cand))]  # ***
    return [tbl for tbl in table_candidates if _no_conflicts(tbl)]


# TODO: Reconcile this version with the similar method in CayleyTable.
def is_table_associative(table):
    """Determine whether the table supports an associative operation."""
    # result = True
    elements = table[0]  # The first row should correspond to the elements of a group
    for a in elements:
        for b in elements:
            for c in elements:
                ab = table[a][b]
                bc = table[b][c]
                if not (table[ab][c] == table[a][bc]):
                    # result = False
                    # break
                    return False
    # return result
    return True


def tables_to_groups(tables, identity_name="e", elem_name="a"):
    """Given a list of multiplication tables, all the same size, turn them into a list of groups."""
    order = len(tables[0])
    groups = []
    for j in range(len(tables)):
        gname = f"G{j}"
        desc = f"Group {j} of order {order}"
        ename = elem_name + str(j) + "_"
        elements = [identity_name + str(j)] + [f"{ename}" + str(i) for i in range(1, order)]
        groups.append(Group(gname, desc, elements, tables[j]))
    return groups


def get_integer_form(elem_list):
    """For an element list like ['e1', 'a1_2', 'a1_1', 'a1_3'],
    return the integer 213, i.e., the 'subscripts' of the elements that
    follow the identity element. Used by get_int_forms."""
    return int(''.join(map(lambda x: x.split("_")[1], elem_list[1:])))


def get_int_forms(ref_group, isomorphisms):
    """Return a list of integer forms ('permutations') for a list of isomorphisms,
    i.e., mappings, based on a reference group."""
    return [get_integer_form([iso[elem] for elem in ref_group.elements])
            for iso in isomorphisms]


