# =========
#   Group
# =========

import numpy as np
from pprint import pprint

from finalg.monoid import Monoid
# from main import make_finite_algebra
from finalg.about_subalgebras import partition_into_isomorphic_lists
from finalg.utilities import yes_or_no


class Group(Monoid):
    """A Group is a Monoid with inverses."""

    def __init__(self, name, description, elements, table, check_inputs=True):
        super().__init__(name, description, elements, table, check_inputs)
        if check_inputs:
            if self.table.has_inverses():
                self._inverses = self.create_inverse_lookup_dict()
            else:
                raise ValueError("CHECK INPUTS: Table has insufficient inverses")
        else:
            self._inverses = self.create_inverse_lookup_dict()

    def __truediv__(self, normal_subgroup):
        """Return the quotient group based on a given normal subgroup.

        Each element of the quotient group will be a representative element from
        one of the subgroup's cosets."""
        return self.quotient_group(normal_subgroup)

    def inverse_mapping(self):
        """Returns a dictionary that maps each element to its inverse."""
        return self._inverses

    def inv(self, element):
        """Return the inverse of an element"""
        return self._inverses[element]

    def sub(self, x, y):
        """Group subtraction:  Return x - y; i.e., x + inv(y)."""
        return self.op(x, self.inv(y))

    def conjugate(self, a, g):
        """Return g * a * inv(g), the conjugate of a with respect to g"""
        return self.op(g, self.op(a, self.inv(g)))

    def commutator(self, a, b):
        """Return [a, b] = a * b * inv(a) * inv(b), the commutator of a & b"""
        return self.op(a, b, self.inv(a), self.inv(b))

    def commutators(self):
        """Return the list of commutators of the group."""
        result = set()
        for a in self:
            for b in self:
                result.add(self.commutator(a, b))
        return result

    def commutator_subalgebra(self):
        """Return the commutator subalgebra (Group, Ring, or Field) of this Group, Ring, or Field."""
        commutators = self.commutators()
        return self.subalgebra_from_elements(self.closure(commutators, include_inverses=True),
                                             name=f"{self.name}_Comm",
                                             desc=f"{self.name} commutator subalgebra")

    def is_normal(self, subgrp):
        """Returns True if the subgroup is normal, otherwise False is returned"""
        for x in self:
            for a in subgrp:
                if not self.conjugate(a, x) in subgrp:
                    return False
        return True

    def trivial_subgroups(self):
        """Return the group's two trivial subgroups."""
        name = f"Subgroup of {self.name}"
        desc = f"Trivial subgroup: {self.description}"
        trivial = Group(name, desc, [self.identity], [[0]])
        return [trivial, self]

    def subgroups(self):
        """Return a list of all subgroups, including trivial subgroups."""
        return self.proper_subalgebras(divisors_only=True, include_inverses=True) + self.trivial_subgroups()

    def unique_proper_subgroups(self, subgroups=None):
        """Return a list of proper subgroups that are unique, up to isomorphism.
        If no subgroups are provided, then they will be derived."""
        if subgroups:
            partitions = partition_into_isomorphic_lists(subgroups)
        else:
            partitions = partition_into_isomorphic_lists(self.proper_subalgebras(divisors_only=True,
                                                                                 include_inverses=True))
        # Return a list of the first subgroups from each sublist of proper subgroups
        return [partition[0] for partition in partitions]

    def quotient_group(self, subgroup):
        """Given a normal subgroup, return the quotient group of this group.
        The elements of the quotient group will be representative elements from
        cosets, prefixed with '~'."""
        from finalg.make_finite_algebra import make_finite_algebra

        def index_of_coset(elem, _cosets):
            """Given an element of an algebra and a list of cosets, find the position of the coset
            that contains the element in the list of cosets."""
            index = None
            for coset in _cosets:
                if elem in coset:
                    index = _cosets.index(coset)
            return index

        if self.is_normal(subgroup):
            cosets = list(self.left_cosets(subgroup))
        else:
            raise ValueError(f"{subgroup.name} is not a normal subgroup of {self.name}")

        # Make a list consisting of one representative element from each coset
        elems = [x[0] for x in cosets]
        n = len(elems)
        table = np.zeros((n, n))
        for b in elems:
            b_index = elems.index(b)
            for a in elems:
                a_index = elems.index(a)
                axb = self.op(a, b)
                axb_index = index_of_coset(axb, cosets)
                table[a_index][b_index] = axb_index

        name = f"{self.name}/{subgroup.name}"
        desc = f"Group {self.name} modulo subgroup {subgroup.name}"
        coset_elems = tuple("~" + elem for elem in elems)
        return make_finite_algebra(name, desc, coset_elems, table)

    # This 'about' method differs from the one in FiniteAlgebra in that it prints out
    # more detailed information about elements.
    # TODO: It would be nice to combine the two someday.
    def about(self, max_size=12, max_gens=2, use_table_names=False, show_tables=True, show_elements=True,
              show_generators=False):
        """Print information about the Group."""
        print(f"\n** {self.__class__.__name__} **")
        print(f"Name: {self.name}")
        print(f"Instance ID: {id(self)}")
        print(f"Description: {self.description}")
        print(f"Order: {self.order}")
        print(f"Identity: {repr(self.identity)}")
        print(f"Commutative? {yes_or_no(self.is_commutative())}")
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
        spc = 7
        if show_elements:
            print("Elements:")
            print("   Index   Name   Inverse  Order")
            for elem in self:
                idx_elem = self.elements.index(elem)
                inv_elem = self.inv(elem)
                ord_elem = self.element_order(elem)
                # repr is used below to make strings explicit by printing out their quotation marks.
                print(f"{idx_elem :>{spc}} {repr(elem) :>{spc}} {repr(inv_elem) :>{spc}} {ord_elem :>{spc}}")
        size = len(self.elements)
        if show_tables:
            if size <= max_size:
                if use_table_names:
                    print(f"Cayley Table (showing names):")
                    pprint(self.table.to_list_with_names(self.elements))
                else:
                    print(f"Cayley Table (showing indices):")
                    pprint(self.table.tolist())
            else:
                print(f"{self.__class__.__name__} order is {size} > {max_size}, so no table is printed.")
        return str(self)

# def left_cosets(group, subgroup):
#     """Returns an iterator that returns lists of left cosets."""
#     return map(lambda s: sorted(list(s)),
#                {frozenset([group.op(x, y) for y in subgroup])
#                 for x in group})
#
# def right_cosets(group, subgroup):
#     """Returns an iterator that returns lists of right cosets."""
#     return map(lambda s: sorted(list(s)),
#                {frozenset([group.op(y, x) for y in subgroup])
#                 for x in group})

