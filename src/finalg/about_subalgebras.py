# ##########################################################################
# The following functions are not Group methods, even though some of them
# have argument names that sound like they're groups, and descriptions that
# refer to groups.
# ##########################################################################

from functools import reduce

# def partition_into_isomorphic_lists(list_of_groups):
#     """Partition the list of groups into sub-lists of groups that are isomorphic to each other.
#     The purpose of this function is to operate on the proper subgroups of a group to determine
#     the unique subgroups, up to isomorphism.
#     """
#     def iso_and_not_iso(gp, gps):
#         """Partition the list of groups, gps, into two lists, those that are isomorphic to gp
#         and those that are not."""
#         iso_to_grp = []
#         not_iso_to_grp = []
#         for g in gps:
#             if gp.fast_isomorphic(g):
#                 iso_to_grp.append(g)
#             else:
#                 not_iso_to_grp.append(g)
#         return iso_to_grp, not_iso_to_grp
#
#     def aux(result, remainder):
#         """Recursively partition 'remainder' into lists that are isomorphic to its first member of the
#         remainder list and those that are not.  Then, put those that are isomorphic to the first member
#         into the 'result' list, and recurse on the remainder.
#         """
#         if len(remainder) == 0:
#             return result
#         else:
#             first = remainder[0]
#             rest = remainder
#             iso_to_first, not_iso_to_first = iso_and_not_iso(first, rest)
#             result.append(iso_to_first)
#             return aux(result, not_iso_to_first)
#
#     return aux([], list_of_groups)

# def partition_into_isomorphic_lists(list_of_groups):
#     """Partition the list of groups into sub-lists of groups that are isomorphic to each other.
#     The purpose of this function is to operate on the proper subgroups of a group to determine
#     the unique subgroups, up to isomorphism.
#
#     Uses 'fast_isomorphic' rather than the brute-force 'isomorphic', since the latter tries
#     every one of other.order! bijections and is intractable once subalgebra orders climb
#     past roughly 8-9 (e.g. the order-12 subgroups of A5).
#     """
#     def iso_and_not_iso(gp, gps):
#         """Partition the list of groups, gps, into two lists, those that are isomorphic to gp
#         and those that are not."""
#         iso_to_grp = []
#         not_iso_to_grp = []
#         for g in gps:
#             if gp.fast_isomorphic(g):
#                 iso_to_grp.append(g)
#             else:
#                 not_iso_to_grp.append(g)
#         return iso_to_grp, not_iso_to_grp
#
#     def aux(result, remainder):
#         """Recursively partition 'remainder' into lists that are isomorphic to its first member
#         and those that are not, then recurse on the remainder."""
#         if len(remainder) == 0:
#             return result
#         else:
#             first = remainder[0]
#             rest = remainder[1:]  # was 'remainder', which wastefully compared 'first' to itself
#             iso_to_first, not_iso_to_first = iso_and_not_iso(first, rest)
#             result.append([first] + iso_to_first)
#             return aux(result, not_iso_to_first)
#
#     return aux([], list_of_groups)


def partition_into_isomorphic_lists(algebras):
    """
    Groups a list of Algebra instances into equivalence classes based on isomorphism,
    and sorts the resulting classes by their Algebra order from smallest to largest.

    Parameters:
    algebras (list): A list of Algebra class instances.

    Returns:
    list: A list of lists, where each inner list represents an equivalence class,
          sorted ascending by the 'order' property of the algebras.
    """
    classes = []

    for current_algebra in algebras:
        placed = False
        # Check if current_algebra belongs to an existing equivalence class
        for eq_class in classes:
            # Check against the first representative of the class
            representative = eq_class[0]
            if current_algebra.fast_isomorphic(representative):
                eq_class.append(current_algebra)
                placed = True
                break

        # If it doesn't match any existing class, it starts its own new class
        if not placed:
            classes.append([current_algebra])

    # Sort the list of classes based on the 'order' of the first item in each class
    classes.sort(key=lambda eq_class: eq_class[0].order)

    return classes


def about_isomorphic_partition(alg, part, identify_small=True):
    """Print a summary of a particular partition of isomorphic subalgebras of an algebra.

    If identify_small is True (the default) and the subalgebras are groups of order 20 or less, then
    the summary also names the small group (see finalg.small_algebras) that they are isomorphic to.
    """
    size = len(part)

    if size == 0:
        raise ValueError("A partition must have at least one member.")

    # All the algebras in a partition are isomorphic to each other,
    # so, get (most) properties from the first algebra in the partition
    sub0 = part[0]
    classname = f"{sub0.__class__.__name__}"
    order = sub0.order

    # If the subalgebras are isomorphic to a small group, then the summary says which one.
    small_note = ""
    if identify_small:
        from finalg.small_algebras import identify_small_group  # Imported here to avoid circular imports
        small = identify_small_group(sub0)
        if small is not None:
            each = "each " if size > 1 else ""
            small_note = f"   ({each}isomorphic to the small group {small.name})"

    comm = "Isomorphic "
    comm1 = ""
    if sub0.is_commutative():
        comm = "Isomorphic Commutative "
        comm1 = "Commutative "

    # See if there are different identity elements for each algebra in the partition
    identities = False
    single_id = False
    if sub0.has_identity():
        identities = {sub.identity for sub in part}
        if len(identities) == 1:
            single_id = sub0.identity

    norm = ""
    if alg.has_inverses() and alg.is_normal(sub0):
        norm = "Normal "

    if size > 1:
        if identities:
            if single_id:
                print(f"{size} {comm}{norm}{classname}s of order {order} with identity '{single_id}':")
                if small_note:
                    print(small_note)
                for sub in part:
                    sub_cname = sub.__class__.__name__
                    print(f"      {sub_cname}: {sub.name}: {sub.elements}")
                print("")
            else:
                print(f"{size} {comm}{norm}{classname}s of order {order}:")
                if small_note:
                    print(small_note)
                for sub in part:
                    sub_cname = sub.__class__.__name__
                    print(f"      {sub_cname}: {sub.name}: {sub.elements} with identity '{sub.identity}'")
                print("")
        else:
            print(f"{size} {comm}{norm}{classname}s of order {order}:")
            if small_note:
                print(small_note)
            for sub in part:
                sub_cname = sub.__class__.__name__
                print(f"      {sub_cname}: {sub.name}: {sub.elements}")
            print("")
    elif size == 1:
        if identities:
            print(f"{size} {comm1}{norm}{classname} of order {order} with identity '{sub0.identity}':")
        else:
            print(f"{size} {comm1}{norm}{classname} of order {order}:")
        if small_note:
            print(small_note)
        sub0_cname = sub0.__class__.__name__
        print(f"      {sub0_cname}: {sub0.name}: {sub0.elements}\n")
    else:
        raise ValueError("A partition must have at least one member.")


def are_n(n):
    """A bit of grammar.  This function returns a string with the appropriate
    singular or plural present indicative form of 'to be', along with 'n'.
    """
    choices = ['are no', 'is 1', f'are {n}']
    if n < 2:
        return choices[n]
    else:
        return choices[2]


def add_s(string, n):
    """Make a string plural by adding an 's' to it, or not, depending on 'n'."""
    if n == 1:
        return string
    else:
        return string + 's'


def about_isomorphic_partitions(alg, partitions, identify_small=True):
    """Print a summary of the isomorphic partitions of an algebra. If identify_small is True (the default),
    then each partition of groups of order 20 or less also names the small group it is isomorphic to."""
    if len(partitions) != 0:
        n_subs = reduce(lambda x, y: x + y, [len(p) for p in partitions])
        n_parts = len(partitions)
        print(f"\nSubalgebras of {alg.name} : {alg.description}\n")
        what = f"  There {are_n(n_parts)} unique proper {add_s('subalgebra', n_parts)}, up to isomorphism, "
        out_of = f"out of {n_subs} total subalgebras."
        print(what + out_of)
        print(f"  as shown below:\n")
        for partition in partitions:
            about_isomorphic_partition(alg, partition, identify_small)
    else:
        print("There are no proper subalgebras.")


def about_subalgebras(alg, identify_small=True):
    """A convenience function that finds and summarizes all proper subalgebras
    of the input FiniteAlgebra.  The list of isomorphic partitions is
    returned and a summary of it is printed out.

    If identify_small is True (the default), then, for each set of isomorphic subgroups of
    order 20 or less, the summary names the small group that they are isomorphic to, e.g., Z4,
    V4, or D3 (see finalg.small_groups). Set identify_small to False to leave this out.
    """
    alg_subs = alg.proper_subalgebras()
    partitions = partition_into_isomorphic_lists(alg_subs)
    about_isomorphic_partitions(alg, partitions, identify_small)
    return partitions

def find_isomorphic_subalgebra(algebra, partitions, verbose=False):
    """Given an algebra and the partitions output by the method
    about_subalgebras that was applied to a different algebra,
    find an algebra in the partitions that is isomorphic to the
    given algebra, if one exists. If one is found, return the
    isomorphism (dict) along with the algebra found in the partition.
    Otherwise, return False."""
    iso_grp = None
    iso = False
    n = algebra.order
    for part in partitions:
        if part[0].order == n :
            if verbose:
                print(f"Checking: {part[0].name}")
            iso = algebra.fast_isomorphic(part[0])
            if iso:
                iso_grp = part[0]
                break
    if iso:
        return iso, iso_grp
    else:
        print("Not found.")
        return False

