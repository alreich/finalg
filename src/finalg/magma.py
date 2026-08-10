# =========
#   Magma
# =========

import itertools as it
import numpy as np
import networkx as nx
from pprint import pprint
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

    def closed_subsets_of_elements(self, divisors_only, include_inverses):
        """Return all unique, closed, proper subsets of the algebra's elements.
        This returns a list of lists. Each list represents the elements of a subalgebra.
        If divisors_only is True, then only subalgebras of orders that are divisors of
        self will be examined."""
        closed = set()  # Build the result as a set of sets to avoid duplicates
        all_elements = self.elements
        n = len(all_elements)
        if divisors_only:
            range_ = divisors(n, non_trivial=True)
        else:
            range_ = range(2, n - 1)
        for i in range_:
            # Look at all combinations of elements: pairs, triples, quadruples, etc.
            for combo in it.combinations(all_elements, i):
                # Freezing is required to add a set to a set
                clo = frozenset(self.closure(list(combo), include_inverses))
                if len(clo) < n:  # Don't include closures consisting of all elements
                    closed.add(clo)
        return list(map(lambda x: list(x), closed))

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
    ) -> None:
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
        plt.show()

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


