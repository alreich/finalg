# =========
#   Magma
# =========

import itertools as it
from collections import Counter
import numpy as np
import networkx as nx
from pprint import pprint
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

from finalg.finite_operator import FiniteOperator
from finalg.finite_algebra_ABC import FiniteAlgebra
# from main import make_finite_algebra
# from ring import Ring
from finalg.utilities import yes_or_no
from finalg.my_math import divisors

class Magma(FiniteAlgebra):
    """A Magma is a finite algebra with a binary operation that returns a unique value, in the algebra,
    for all pairs in the cross-product of the algebra's set of elements with itself.  With a binary
    operation, we can compute the direct product of two or more algebras.  Also, we can check to see
    if two Magmas are isomorphic."""

    def __init__(self, name, description, elements, table):
        super().__init__(name, description, elements, table)
        self.op = FiniteOperator(self.elements, self.identity, self.table)
        self._dp_delimiter = ':'  # name delimiter used when creating Direct Products

    def _key(self):
        return tuple([self.elements, tuple(tuple(row) for row in self.table.tolist())])

    def __hash__(self):
        return hash(self._key())

    def __eq__(self, other):
        if isinstance(other, Magma):
            return self._key() == other._key()
        else:
            return NotImplemented

    def direct_product_delimiter(self, delimiter=None):
        """If no input, then the current direct product element name delimiter will be returned (default is ':').
        Otherwise, if a string is input (e.g., "-") it will become the new delimiter for direct product element
        names, and then it will be returned."""
        if delimiter:
            self._dp_delimiter = delimiter
            return self._dp_delimiter
        else:
            return self._dp_delimiter

    def __mul__(self, other):  # Direct Product of two Magmas
        """Return direct product of this algebra with the `other` algebra."""
        from finalg.make_finite_algebra import make_finite_algebra
        dp_name = f"{self.name}_x_{other.name}"
        dp_description = "Direct product of " + self.name + " & " + other.name
        dp_element_names = list(it.product(self.elements, other.elements))  # Cross product
        dp_mult_table = list()
        for a in dp_element_names:
            dp_mult_table_row = list()  # Start a new row
            for b in dp_element_names:
                dp_mult_table_row.append(dp_element_names.index((self.op(a[0], b[0]), other.op(a[1], b[1]))))
            dp_mult_table.append(dp_mult_table_row)  # Append the new row to the table
        return make_finite_algebra(dp_name,
                                   dp_description,
                                   list([f"{elem[0]}{self._dp_delimiter}{elem[1]}" for elem in dp_element_names]),
                                   dp_mult_table)

    # TODO: Implement Semi-Direct Product

    def __pow__(self, n, modulo=None):
        """Return the direct product of this algebra with itself, n times.
        Ignore the modulo argument for algebras."""
        result = self
        if isinstance(n, int) and n > 0:
            for _ in range(n - 1):
                result = result * self
        else:
            raise ValueError(f"Error: n = {n}. Power must be a positive integer.")
        return result

    def element_to_power(self, elem, n, left_associative=True):
        """Return the n_th power of the given element. n must be an integer. If n == 0 and an
        identity element exists, then it will be returned; otherwise, a ValueError is raised.
        If n < 0, and the algebra has inverses, then the inverse of the element raised to the
        absolute value of the power is returned, e.g., b^-4 = inv(b^4). If n < 0 and the algebra
        does not have inverses, then a ValueError is raised. For non-associative algebras (Magmas),
        the default is for products to be associated from the left, e.g., b^4 = ((b * b) * b) * b.
        Set left_associative to False, to associate from the right, instead."""
        result = elem
        if elem in self.elements:
            if isinstance(n, int):  # Or, raise an exception if not an integer
                if n == 0:
                    if self.has_identity():  # Or, raise an exception if no identity exists
                        result = self.identity
                    else:
                        raise ValueError(f"Error: n = {n}. {self.name} does not have an identity element.")
                elif n > 0:
                    for _ in range(n - 1):
                        if left_associative:
                            result = self.op(result, elem)
                        else:
                            result = self.op(elem, result)
                elif n < 0:  # Or, raise an exception if inverses do not exist
                    if self.has_inverses():
                        result = self.inv(self.element_to_power(elem, abs(n), left_associative))
                    else:
                        raise ValueError(f"Error: n = {n}. But, {self.name} does not have inverses.")
            else:
                raise ValueError(f"Error: n = {n}. The power must be an integer.")
        else:
            raise ValueError(f"Error: {elem} is not an element of {self.name}.")
        return result

    def reorder_elements(self, reordered_elements):
        """Return a new group made from this one with the elements reordered."""
        n = self.order
        if n == len(reordered_elements):
            new_table = np.full((n, n), 0)
            for row in range(n):
                for col in range(n):
                    prod = self.op(reordered_elements[row], reordered_elements[col])
                    new_table[row, col] = reordered_elements.index(prod)
            new_name = str(self.name) + '_REORDERED'
            new_desc = str(self.description) + ' (elements reordered)'
            return self.__class__(new_name, new_desc, tuple(reordered_elements), new_table)
        else:
            raise ValueError(f"There are {len(reordered_elements)} reordered elements.  There should be {n}.")

    def make_element_mappings(self, other):
        """Returns a list of mappings (dictionaries) of this algebra's elements to all
        possible permutations of other's elements.  The orders of self and other must
        be equal."""
        return [dict(zip(self.elements, perm)) for perm in it.permutations(other.elements)]

    def is_isomorphic_mapping(self, other, mapping):
        """Returns True if the input mapping from this algebra to the other algebra is isomorphic."""
        return all([mapping[self.op(x, y)] == other.op(mapping[x], mapping[y])
                    for x in self.elements for y in self.elements])

    def isomorphic(self, other):
        """If there is a mapping from elements of this algebra to the other algebra's elements,
        return it; otherwise return False."""
        if (self.__class__.__name__ == other.__class__.__name__) and (self.order == other.order):
            maps = self.make_element_mappings(other)
            for mp in maps:
                if self.is_isomorphic_mapping(other, mp):
                    return mp
            return False
        else:
            return False

    # ------------------------------------------------------------------
    # Faster isomorphism testing via generating sets ("subalgebras")
    # ------------------------------------------------------------------
    #
    # 'isomorphic', above, is brute force: it tries every one of
    # other.order! bijections between the two element sets. That is fine
    # for small algebras, but it is hopeless once the order climbs past
    # roughly 8 or 9.
    #
    # 'fast_isomorphic' instead exploits the fact that an isomorphism is
    # completely determined by where it sends a *generating set* of this
    # algebra: once the images of the generators are chosen, the image of
    # every other element is forced, because the generators' closure
    # (i.e. the chain of subalgebras built up by repeatedly multiplying
    # known elements together) is the whole algebra. So instead of
    # searching over all bijections of the full element set, we only
    # search over bijections of a (usually much smaller) generating set,
    # and "grow" each candidate the same way 'closure' grows a
    # subalgebra -- except here we grow the mapping alongside it,
    # checking consistency at every step and bailing out of bad
    # candidates as early as possible. Candidate images are also
    # restricted, up front, to elements that share the same invariants
    # (row/column "fingerprint", idempotency, element order, etc.),
    # which are necessarily preserved by any isomorphism.

    def _row_profile(self, x):
        """Return a sorted tuple describing, for element x, the sizes of the
        'fibers' of the map y -> self.op(x, y); i.e., for each result value,
        how many y produce it. This profile (ignoring which result value goes
        with which count) is preserved by any isomorphism: if f is an
        isomorphism, then |{y : x*y = z}| == |{y' : f(x)*y' = f(z)}| for all z."""
        counts = Counter(self.op(x, y) for y in self.elements)
        return tuple(sorted(counts.values()))

    def _col_profile(self, x):
        """Column analog of _row_profile: the fiber-size profile of y -> self.op(y, x)."""
        counts = Counter(self.op(y, x) for y in self.elements)
        return tuple(sorted(counts.values()))

    def _element_invariant(self, x):
        """Return a tuple of properties of element x that must be shared by its
        image under any isomorphism. Used to prune candidate images before
        searching for an isomorphism, rather than trying every element."""
        inv = [self._row_profile(x), self._col_profile(x), self.op(x, x) == x]
        if hasattr(self, 'element_order'):  # Monoids (& subclasses) can compute this cheaply
            inv.append(self.element_order(x))
        return tuple(inv)

    def _element_invariants(self):
        """Return a dict mapping each element of this algebra to its _element_invariant."""
        return {x: self._element_invariant(x) for x in self.elements}

    def _smallest_generating_set(self):
        """Return a list of elements of this algebra whose closure, under the
        algebra's operation, is the entire algebra, preferring the smallest
        such set available (as found by 'generators')."""
        gens = self.generators()
        if len(gens) == 0:
            return list(self.elements)  # Degenerate fallback; shouldn't normally happen.
        first = gens[0]
        # When the algebra is cyclic, 'generators' returns a flat list of individual
        # elements (each one generates the algebra alone); otherwise, it returns a
        # list of tuples, each of which is a minimal generating set.
        if isinstance(first, (tuple, list)):
            return list(first)
        else:
            return [first]

    def _extend_mapping(self, other, partial_map):
        """Given a partial mapping (dict) whose keys generate this algebra (i.e.,
        repeatedly applying self.op to known key/value pairs eventually reaches
        every element), grow it into a mapping over all of self's elements by
        propagating products: whenever x1 -> y1 and x2 -> y2 are known, x1*x2 must
        map to y1*y2. If this ever conflicts with an existing entry, or the result
        isn't a bijection covering every element, return None. Otherwise, return
        the completed mapping if it verifies as an isomorphism, else None."""
        mapping = dict(partial_map)
        n = self.order
        growing = True
        while growing and len(mapping) < n:
            growing = False
            items = list(mapping.items())
            for x1, y1 in items:
                for x2, y2 in items:
                    x3 = self.op(x1, x2)
                    y3 = other.op(y1, y2)
                    if x3 in mapping:
                        if mapping[x3] != y3:
                            return None  # Inconsistent: not a valid candidate.
                    else:
                        mapping[x3] = y3
                        growing = True
        if len(mapping) != n or len(set(mapping.values())) != n:
            return None  # Generators didn't reach every element, or map isn't a bijection.
        return mapping if self.is_isomorphic_mapping(other, mapping) else None

    def fast_isomorphic(self, other, verbose=False):
        """A faster alternative to 'isomorphic' for determining whether this algebra
        and 'other' are isomorphic, and if so, returning the mapping (a dict) between
        their elements.

        Instead of brute-force trying all other.order! bijections, this method:
          1. Rules out obviously non-isomorphic algebras cheaply (class, order,
             whether an identity exists, commutativity).
          2. Computes an isomorphism-invariant signature for every element of both
             algebras (see '_element_invariant'), and bails out immediately if the
             multiset of signatures doesn't match between the two algebras.
          3. Finds a small generating set for this algebra -- a set of elements
             whose repeated products ("closure") reach every element, i.e., build up
             a chain of subalgebras that ends at the whole algebra.
          4. Tries only the images for those generators that share the right
             signature, and, for each candidate, extends it to a full mapping by
             propagating products the same way 'closure' grows a subalgebra,
             checking consistency at each step and abandoning bad candidates as
             soon as a conflict appears.

        Because step 3 usually finds a generating set that is far smaller than the
        whole algebra, and step 2 shrinks the candidate images for each generator,
        the number of candidates actually tried is normally a tiny fraction of
        other.order!, and inconsistent candidates tend to fail fast rather than
        being checked against the whole table. Set verbose=True to see a little
        information about the search as it runs.

        Returns the mapping (dict) if the algebras are isomorphic; otherwise False.
        """
        if (self.__class__.__name__ != other.__class__.__name__) or (self.order != other.order):
            return False
        if self.has_identity() != other.has_identity():
            return False
        if self.is_commutative() != other.is_commutative():
            return False

        n = self.order
        if n == 1:
            return {self.elements[0]: other.elements[0]}

        inv0 = self._element_invariants()
        inv1 = other._element_invariants()

        if Counter(inv0.values()) != Counter(inv1.values()):
            if verbose:
                print("Element invariants don't match between the two algebras; not isomorphic.")
            return False

        if self.has_identity() and inv0[self.identity] != inv1[other.identity]:
            return False

        gens = self._smallest_generating_set()

        # Group other's elements by invariant, so candidate images can be looked up directly.
        by_invariant = {}
        for e in other.elements:
            by_invariant.setdefault(inv1[e], []).append(e)
        candidate_lists = [by_invariant.get(inv0[g], []) for g in gens]

        if any(len(cands) == 0 for cands in candidate_lists):
            return False  # A generator has no element in 'other' with a matching invariant.

        if verbose:
            print(f"Generating set (size {len(gens)}): {gens}")
            print(f"Candidate image counts per generator: {[len(c) for c in candidate_lists]}")

        tried = 0
        for images in it.product(*candidate_lists):
            if len(set(images)) != len(images):
                continue  # Images of distinct generators must be distinct.
            tried += 1
            partial = dict(zip(gens, images))
            if self.has_identity():
                partial[self.identity] = other.identity
            mapping = self._extend_mapping(other, partial)
            if mapping is not None:
                if verbose:
                    print(f"Isomorphism found after trying {tried} generator-image assignment(s).")
                return mapping

        if verbose:
            print(f"Tried {tried} generator-image assignment(s); no isomorphism found.")
        return False

    def closure(self, subset_of_elements, include_inverses):
        """Given a subset (in list form) of the group's elements (name strings),
        return the smallest possible set of elements containing the subset
        that is closed under the algebra's operation(s).  If include_inverses
        is True and the algebra has inverses, then they will be added to the
        closure."""

        result = set(subset_of_elements)

        # Include inverses, maybe.
        if include_inverses and self.has_inverses():
            for elem in subset_of_elements:
                result.add(self.inv(elem))

        # Add the products (sums, if rings) of all possible pairs
        for pair in it.product(result, result):
            result.add(self.op(*pair))

        # For rings, add the products of all possible pairs
        # if isinstance(self, Ring):
        if hasattr(self, 'mult'):
            for pair in it.product(result, result):
                result.add(self.mult(*pair))

        # If the input set of elements increased, recurse ...
        if len(result) > len(subset_of_elements):
            return self.closure(list(result), include_inverses)

        # ...otherwise, stop and return the result
        else:
            return list(result)

    # def closed_subsets_of_elements(self, divisors_only, include_inverses):
    #     """Return all unique, closed, proper subsets of the algebra's elements.
    #     This returns a list of lists. Each list represents the elements of a subalgebra.
    #     If divisors_only is True, then only subalgebras of orders that are divisors of
    #     self will be examined."""
    #     closed = set()  # Build the result as a set of sets to avoid duplicates
    #     all_elements = self.elements
    #     n = len(all_elements)
    #     if divisors_only:
    #         range_ = divisors(n, non_trivial=True)
    #     else:
    #         range_ = range(2, n - 1)
    #     for i in range_:
    #         # Look at all combinations of elements: pairs, triples, quadruples, etc.
    #         for combo in it.combinations(all_elements, i):
    #             # Freezing is required to add a set to a set
    #             clo = frozenset(self.closure(list(combo), include_inverses))
    #             if len(clo) < n:  # Don't include closures consisting of all elements
    #                 closed.add(clo)
    #     return list(map(lambda x: list(x), closed))

    def closed_subsets_of_elements(self, divisors_only=True, include_inverses=True):
        """Return all unique, closed, proper subsets of the algebra's elements.
        This returns a list of lists. Each list represents the elements of a subalgebra.
        If divisors_only is True, only subalgebras of orders that are divisors of self's
        order are included in the result.

        Instead of testing every combination of elements of every candidate size (which is
        combinatorially infeasible -- e.g. C(60, 30) ~ 1.2*10^17 for an order-60 algebra
        like A5), this builds up the lattice of closed subsets from the bottom: starting
        from the closure of each individual element, it repeatedly extends every closed
        set found so far by one more element and re-closes, discovering every closed proper
        subset that way. This works because any closed subset H containing elements
        g1, ..., gk is reachable by adding those elements to the empty set one at a time and
        closing at each step -- closure({g1}) subseteq closure({g1,g2}) subseteq ... subseteq H
        -- and every intermediate closure along the way is itself closed and contained in H,
        so it is one of the sets this search discovers on the way to finding H.

        Reference: this is the "cyclic extension" method for building a subgroup lattice;
        see J. Neubuser, "Untersuchungen des Untergruppenverbandes endlicher Gruppen auf
        einer programmgesteuerten Rechenanlage," Numerische Mathematik 2 (1960): 280-292,
        and D. Holt, B. Eick, E. O'Brien, Handbook of Computational Group Theory, Chapman &
        Hall/CRC, 2005 (subgroup lattice / cyclic extension algorithm).
        """
        all_elements = self.elements
        n = len(all_elements)

        known = set()
        queue = []
        for x in all_elements:
            clo = frozenset(self.closure([x], include_inverses))
            if len(clo) < n and clo not in known:
                known.add(clo)
                queue.append(clo)

        while queue:
            subset = queue.pop()
            for g in all_elements:
                if g in subset:
                    continue
                extended = frozenset(self.closure(list(subset) + [g], include_inverses))
                if len(extended) < n and extended not in known:
                    known.add(extended)
                    queue.append(extended)

        if divisors_only:
            allowed_sizes = set(divisors(n, non_trivial=True))
            known = {s for s in known if len(s) in allowed_sizes}

        return list(map(list, known))

    def subalgebra_from_elements(self, closed_subset_of_elements, name="No name", desc="No description"):
        """Return the algebra constructed from the given closed subset of elements."""
        from finalg.make_finite_algebra import make_finite_algebra
        # Make sure the elements are sorted according to their order in the parent Group (self)
        elements_sorted = sorted(closed_subset_of_elements, key=lambda x: self.elements.index(x))
        table = []
        for a in elements_sorted:
            row = []
            for b in elements_sorted:
                # The table entry is the index of the product in the sorted elements list
                row.append(elements_sorted.index(self.op(a, b)))
            table.append(row)
        # if isinstance(self, Ring):
        if hasattr(self, 'mult'):
            table2 = []
            for c in elements_sorted:
                row2 = []
                for d in elements_sorted:
                    # The table entry is the index of the product in the sorted elements list
                    row2.append(elements_sorted.index(self.mult(c, d)))
                table2.append(row2)
            return make_finite_algebra(name, desc, elements_sorted, table, table2)
        else:
            return make_finite_algebra(name, desc, elements_sorted, table)

    def proper_subalgebras(self, divisors_only=True, include_inverses=True):
        """Return a list of proper subalgebras of the algebra."""
        desc = f"Subalgebra of: {self.description}"
        count = 0
        list_of_subalgebras = []
        for closed_element_set in self.closed_subsets_of_elements(divisors_only, include_inverses):
            name = f"{self.name}_subalgebra_{count}"
            count += 1
            list_of_subalgebras.append(self.subalgebra_from_elements(closed_element_set, name, desc))
        return list_of_subalgebras

    def generates(self, set_of_elems):
        """Returns True if a set of one or more elements generates the algebra,
        otherwise False is returned.
        """
        clo = self.closure(set_of_elems, include_inverses=False)
        return set(clo) == set(self.elements)

    # TODO: Review this method and the is_cyclic method for issues (see notebook)
    def generators(self, start_of_range=1):
        """If the algebra is cyclic, then a list of individual elements that each
        generate the algebra is returned; otherwise, a list of lists of elements,
        is returned, where each sublist generates the algebra. This method looks
        for the smallest sets of elements that can generate the group. It stops
        looking once it finds all small sets of elements of a given size.
        """
        gens = list()
        n = None  # define n outside the scope of the loop below
        for k in range(start_of_range, self.order + 1):
            n = k  # Save for test later
            combos = list(it.combinations(self.elements, k))
            stop = False
            for combo in combos:
                if self.generates(combo):
                    gens.append(combo)
                    stop = True
            if stop:
                break
        if n == 1:
            return [gen[0] for gen in gens]  # The grp is cyclic
        else:
            return gens  # The grp is not cyclic

    def get_single_generator_set(self):
        """A convenience function that returns the first list of generators
        from the list of all generators returned by the method 'generators'.
        This function also makes sure that the generator set returned is a
        list, and not just a single generator element."""
        gens = self.generators()
        gen = gens[0]
        if isinstance(gen, tuple) or isinstance(gen, tuple):
            return list(gen)
        else:
            return list(gen)

    def is_cyclic(self):
        """Returns False if this algebra is not cyclic; otherwise a list of elements
        is returned, where each one can generate the entire algebra."""
        elemset = set(self.elements)
        gens = [x for x in elemset if set(self.closure([x], False)) == elemset]
        if len(gens) == 0:
            return False
        else:
            return gens

    def center(self):
        """Return the list of elements that commute with every element of the algebra.
        In Pinter's book, chapter 5, exercise D3, the 'center' is defined for Groups,
        but the definition also works for any Magma."""
        return [a for a in self if all([self.op(a, x) == self.op(x, a) for x in self])]

    def center_algebra(self, verbose=False):
        """Return the subalgebra that is the center of this algebra.  If the center is part of a
        Semigroup, then (due to associativity) it will be closed wrt the Semigroup operation,
        and hence form a sub-semigroup, but the center of a Magma will not necessarily be closed.
        Note also that, if the algebra is commutative, then the entire algebra is its center."""
        ctr = self.center()

        if len(ctr) == 0:  # If there is no center...
            if verbose:
                print(f"{self} does not have a Center.")
            return None

        # Being lazy (or careful) and checking closure, regardless of the type of algebra
        elif set(ctr) == set(self.closure(ctr, False)):
            return self.subalgebra_from_elements(ctr, self.name + '_CENTER', 'Center of ' + self.name)
        else:
            if verbose:
                print(f"The Center of {self} is not closed.")
            return None

    # def is_division_algebra(self, verbose=False):
    #     """Return True if, for every a & b in the algebra, there are unique x and y in the algebra
    #     such that ax=b and ya=b. Otherwise, return False. Set verbose to True to see intermediate
    #     calculations."""
    #     if verbose:
    #         print(f"\n{self}\n")
    #     result = True
    #     elems = self.elements
    #     n_sqr = self.order ** 2
    #     count = 0  # number of successes
    #     for ab in it.product(elems, elems):
    #         a = ab[0]
    #         b = ab[1]
    #         ab_ok = False
    #         for xy in it.product(elems, elems):
    #             x = xy[0]
    #             y = xy[1]
    #             if self.op(a, x) == b and self.op(y, a) == b:
    #                 count += 1
    #                 if verbose:
    #                     print(f"{ab} & {xy}")
    #                 ab_ok = True
    #                 # break
    #         if not ab_ok:
    #             result = False
    #             if verbose:
    #                 print(f"{ab} fail")
    #     if verbose:
    #         print(f"Number of successes, {count}, should equal {n_sqr}")
    #     if result:
    #         if count == n_sqr:
    #             return True
    #         elif count > n_sqr:
    #             if verbose:
    #                 print(f"Count of {count} > {n_sqr} means some cancellations are not unique.")
    #             return False
    #         else:
    #             raise Exception(f"A True result with count {count} < {n_sqr} means something went wrong.")
    #     else:
    #         return False

    def _element_pairs_where_table_equals(self, cayley_table, elem_name):
        """Utility function that returns all pairs of elements where the cayley_table entries
         are equal to elem."""
        elems = self.elements
        index = elems.index(elem_name)
        pairs = cayley_table.table_entries_where_equal_to(index)
        return [(elems[pair[0]], elems[pair[1]]) for pair in pairs]

    def element_pairs_where_sum_equals(self, elem_name):
        """Return all pairs of elements where the sums are equal to elem_name. The 'sum' here
        refers to the binary operation of a Magma, Semigroup, Group, or to the additive binary
        operation of a Ring or Field.
        """
        return self._element_pairs_where_table_equals(self.table, elem_name)

    def left_cosets(self, subalgebra):
        """Returns an iterator that returns lists of left cosets."""
        return map(lambda s: sorted(list(s)),
                   {frozenset([self.op(x, y) for y in subalgebra])
                    for x in self})

    def right_cosets(self, subalgebra):
        """Returns an iterator that returns lists of right cosets."""
        return map(lambda s: sorted(list(s)),
                   {frozenset([self.op(y, x) for y in subalgebra])
                    for x in self})

    def make_cayley_graph(self, generators) -> nx.MultiDiGraph:
        """
        Build a directed multigraph whose nodes are group elements and
        whose edges are labeled by generators. The Edge attribute, gen_idx,
        stores the index, in the list of generators, of the generator
        represented by the edge.
        """
        # Give every element a canonical integer label
        ilabel = {elem: i for i, elem in enumerate(self.elements)}

        # Create the graph and add its nodes
        cayley_graph = nx.MultiDiGraph()
        cayley_graph.add_nodes_from(range(self.order))
        nx.set_node_attributes(cayley_graph,
                               {i: str(elem) for i, elem in enumerate(self.elements)},
                               "label")

        # Create graph edges & label them according to the generator each one represents
        for elem in self.elements:
            for i, gen in enumerate(generators):
                product = self.op(gen, elem)
                cayley_graph.add_edge(ilabel[elem], ilabel[product], gen_idx=i)

        return cayley_graph

    def draw_cayley_diagram(
            self,
            generators=None,
            layout: str = "spring",  # "spring" | "circular" | "shell" | "spectral"
            legend_loc: str = "upper left",
            node_size: int = 800,
            font_size: int = 12,
            figsize: tuple[int, int] = (9, 7),
            show: bool = True,
    ) -> plt.Figure:
        """
        Draw the Cayley diagram of `group` with respect to `generators`.

        Parameters
        ----------
        self      : list of Permutation objects that form a finite group.
        generators : explicit generating set; auto-detected if None.
        layout     : NetworkX layout algorithm to use.
        node_size  : matplotlib node size.
        font_size  : label font size.
        figsize    : figure size in inches.
        legend_loc : location of legend in drawing
        show       : if True (default), display the figure via `plt.show()`,
                     but only when running under an interactive Matplotlib
                     backend. Under a non-interactive backend (e.g. "Agg",
                     as used by automated tests), `plt.show()` is skipped to
                     avoid the "FigureCanvasAgg is non-interactive" UserWarning.
                     Set to False to always skip displaying the figure.

        Returns
        -------
        matplotlib.figure.Figure
            The created figure, so callers (including tests) can inspect it
            without needing the figure to be shown.
        """

        # Color palette for up to 8 generators
        _PALETTE = [
            "#e05252",  # red
            "#4a90d9",  # blue
            "#5cb85c",  # green
            "#f0ad4e",  # orange
            "#9b59b6",  # purple
            "#1abc9c",  # teal
            "#e74c3c",  # crimson
            "#34495e",  # dark-slate
        ]

        if generators is None:
            generators = self.get_single_generator_set()

        G = self.make_cayley_graph(generators)
        n_gen = len(generators)
        colors = _PALETTE[:n_gen]

        # ── layout ────────────────────────────────────────────────────────────────
        layouts = {
            "spring": lambda: nx.spring_layout(G, seed=42, k=2.5 / self.order ** 0.5),
            "circular": lambda: nx.circular_layout(G),
            "shell": lambda: nx.shell_layout(G),
            "spectral": lambda: nx.spectral_layout(G),
        }
        pos = layouts.get(layout, layouts["spring"])()

        # ── draw ──────────────────────────────────────────────────────────────────
        fig, ax = plt.subplots(figsize=figsize)
        title = self.name + " : " + self.description
        ax.set_title(title, fontsize=14, fontweight="bold", pad=16)
        ax.axis("off")

        # Nodes
        nx.draw_networkx_nodes(
            G, pos, ax=ax,
            node_size=node_size,
            node_color="#f8f8f2",
            edgecolors="#333333",
            linewidths=1.5,
        )
        # Labels
        labels = nx.get_node_attributes(G, "label")
        nx.draw_networkx_labels(G, pos, labels, ax=ax, font_size=font_size, font_color="#222222")

        # Edges — one pass per generator so we can color them
        for gen_idx in range(n_gen):
            edge_list = [
                (u, v)
                for u, v, d in G.edges(data=True)
                if d["gen_idx"] == gen_idx
            ]
            nx.draw_networkx_edges(
                G, pos, ax=ax,
                edgelist=edge_list,
                edge_color=colors[gen_idx],
                arrows=True,
                arrowsize=18,
                arrowstyle="-|>",
                connectionstyle="arc3,rad=0.12",
                width=2.0,
                min_source_margin=18,
                min_target_margin=18,
            )
        # Legend
        handles = [
            mpatches.Patch(color=colors[i], label=f"× {generators[i]}")
            for i in range(n_gen)
        ]
        ax.legend(
            handles=handles,
            title="Generators",
            loc=legend_loc,
            fontsize=9,
            framealpha=0.85,
        )
        plt.tight_layout()

        # Only pop up the figure under an interactive backend, and only if the
        # caller hasn't opted out. Under a non-interactive backend (e.g. "Agg",
        # as set by the test suite), plt.show() would just emit
        # "FigureCanvasAgg is non-interactive, and thus cannot be shown", so
        # skip it there.
        if show and matplotlib.get_backend().lower() != "agg":
            plt.show()

        return fig

    # This 'about' method differs from the one in Groups in that it does not print out
    # as much detailed information about elements.
    # TODO: Combine the 'about' method, below, with the one in Groups.
    def about(self, max_size=12, max_gens=2, use_table_names=False, show_tables=True,
              show_elements=True, show_generators=False):
        """Prints out information about the algebra. Tables larger than
        max_size are not printed out."""
        print(f"\n** {self.__class__.__name__} **")
        print(f"Name: {self.name}")
        print(f"Instance ID: {id(self)}")
        print(f"Description: {self.description}")
        print(f"Order: {self.order}")
        if self.identity is None:
            print("Identity: None")
        else:
            print(f"Identity: {self.identity}")
        print(f"Associative? {yes_or_no(self.is_associative())}")
        print(f"Commutative? {yes_or_no(self.is_commutative())}")
        # is_cyclic, gens = self.generators()
        single_gens = self.is_cyclic()
        if single_gens:
            print("Cyclic?: Yes")
            print(f"Generators: {single_gens}")
        else:
            print("Cyclic?: No")
            if show_generators:
                gens = self.generators(2)
                num_gens = len(gens)
                gens_sorted = sorted(gens)
                if num_gens > max_gens:
                    print(f"Generators: {gens_sorted[:max_gens]}, plus {num_gens - max_gens} more.")
                elif num_gens == 0:
                    print("Generators: None")
                else:
                    print(f"Generators: {gens_sorted}")
        if show_elements:
            print(f"Elements: {self.elements}")
        print(f"Has Cancellation? {yes_or_no(self.has_cancellation())}")
        print(f"Has Inverses? {yes_or_no(self.has_inverses())}")
        size = len(self.elements)
        if show_tables:
            if size <= max_size:  # Don't print table if too large
                if use_table_names:
                    print(f"Cayley Table (showing names):")
                    pprint(self.table.to_list_with_names(self.elements))
                else:
                    print(f"Cayley Table (showing indices):")
                    pprint(self.table.tolist())
            else:
                print(f"{self.__class__.__name__} order is {size} > {max_size}, so the table is not output.")
        return None


