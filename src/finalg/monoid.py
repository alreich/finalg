# ==========
#   Monoid
# ==========

import numpy as np
import scipy.sparse as sp
import itertools as it

# from ring import Ring
from finalg.semigroup import Semigroup

class Monoid(Semigroup):
    """A Monoid is a Semigroup with an identity element.  With an identity element,
    we can compute element orders.  So, that happens here."""

    def __init__(self, name, description, elements, table, check_inputs=True):
        super().__init__(name, description, elements, table, check_inputs)
        self._element_orders = {elem: 0 for elem in self.elements}  # Cached on first access
        if check_inputs:
            if self.identity is None:
                raise ValueError("CHECK INPUTS: A monoid must have an identity element")

    def element_order(self, element) -> int:
        """Returns the order of the given element within the algebra."""

        def order_aux(elem, prod, order):
            if prod == self.identity:
                return order
            else:
                return order_aux(elem, self.op(prod, elem), order + 1)

        if self._element_orders[element] != 0:  # if already cached, return value
            return self._element_orders[element]
        else:  # else if not cached, compute it, cache it, and return value
            self._element_orders[element] = order_aux(element, element, 1)
            return self._element_orders[element]

    def units(self, return_names=True):
        """Return a sorted list of the Monoid's units.
        By default, the names of elements are returned.
        Setting 'return_names' to False will return element indices instead.
        NOTE: This method is used to compute the units of a Ring.
        """
        # Find the xy-pairs whose product is the Monoid's identity element
        xs, ys = np.where(self.table.table == self.elements.index(self.identity))
        xy_pairs = list(zip(xs, ys))  # e.g., [(1, 2), (2, 1), (5,3), (7, 4), (4, 7)]

        # Collect all x for (x,y) in xy_pairs, if (y,x) is also in xy_pairs
        # e.g., [(1, 2), (2, 1), (5,3), (7, 4), (4, 7)] ==> [1, 2, 4, 7]
        unit_indices = sorted(list({xy[0] for xy in xy_pairs if (xy[1], xy[0]) in xy_pairs}))
        if return_names:
            return [self.elements[index] for index in unit_indices]
        else:
            return unit_indices

    def units_subgroup(self):
        """Return the Unit Subgroup of this algebra.  Makes sense for Monoids or Rings, where
        the multiplicative portion of the Ring is a Monoid.  It will also work for Groups and
        Fields, but will return the entire Group or the entire multiplicative Group of a Field.
        """
        nm = f"{self.name}_Units"
        description = f"Unit subgroup of {self.__class__.__name__}: {self.name}"

        monoid = self
        # if isinstance(self, Ring):
        if hasattr(self, "extract_multiplicative_algebra"):
            monoid = self.extract_multiplicative_algebra()

        return monoid.subalgebra_from_elements(monoid.units(), name=nm, desc=description)

    def regular_representation(self, sparse=""):
        """Given a group, this function returns four things: (1) A dictionary that maps each group
        element to its corresponding regular representation, (2) A (reverse) dictionary that maps each
        regular representation (in the form of a tuple of tuples) back to its corresponding group element,
        (3) A function that maps a group element to its corresponding regular representation matrix, and
        (4) Another function that maps in the opposite direction, from regular representation matrix to
        group element. By default, the matrices are dense arrays. SciPy sparse arrays can be output instead,
        by setting the input variable, "sparse", to one of the following seven strings: "BSR", "COO", "CSC",
        "CSR", "DIA", "DOK", or "LIL". Each of the seven strings corresponds to one of the seven classes of
        sparse array supported by SciPy.
        """
        A = self.elements
        N = self.order

        # Create a list of N Nx1 orthogonal unit vectors
        ident = np.eye(N, dtype=int)  # The NxN identity matrix
        # noinspection PyPep8Naming
        B = [ident[:, [i]] for i in range(N)]  # A list of columns extracted from the identity matrix

        # Create a dictionary that maps group elements to the column vectors created above.
        mapping = dict(zip(A, B))

        # Create a function that takes a group element and returns the corresponding Nx1
        # orthogonal unit vector.
        def V(elem):
            return mapping[elem]

        # Turn a column vector into a tuple, for use as a dict key
        def vector_to_tuple(vec):
            return tuple(map(lambda x: x[0], list(vec)))

        # map vectors to group elements
        inv_mapping = {vector_to_tuple(val): key for key, val in mapping.items()}

        # Given one of the Nx1 orthogonal unit vectors, return the corresponding group element
        def Vinv(vec):
            return inv_mapping[vector_to_tuple(vec)]

        # The seven SciPy sparse array class constructors
        sparse_array_classes = {
            "BSR": sp.bsr_array,
            "COO": sp.coo_array,
            "CSC": sp.csc_array,
            "CSR": sp.csr_array,
            "DIA": sp.dia_array,
            "DOK": sp.dok_array,
            "LIL": sp.lil_array}

        # Create a dictionary that maps each group element to its corresponding
        # regular representation (NxN) matrix.
        reg_rep = dict()
        for k in range(N):
            c_k = np.zeros((N, N))
            for i in range(N):
                for j in range(N):
                    val = np.dot(B[i].transpose(), V(self.op(A[k], Vinv(B[j]))))
                    # c_k[i][j] = np.dot(B[i].transpose(), V(self.op(A[k], Vinv(B[j]))))
                    c_k[i][j] = val.item()
            if sparse in sparse_array_classes:
                reg_rep[A[k]] = sparse_array_classes[sparse](c_k, dtype=int)
            else:
                reg_rep[A[k]] = c_k

        # Create a function that takes a group element and returns the corresponding regular
        # representation matrix, using the dictionary created above.
        def element_to_array(elem):
            return reg_rep[elem]

        # Create a function that turns a 2-dimensional nd.array into a tuple of tuples,
        # for use as a dictionary key. This works for both dense and sparse arrays.
        def array_to_tuple(arr):
            rows, cols = arr.nonzero()
            # return tuple(zip(*arr.nonzero()))
            return tuple(zip(rows, cols))

        # Create a reverse dictionary that maps each regular representation matrix (in tuple form)
        # to its corresponding group element.
        inv_reg_rep = {array_to_tuple(arr): key for key, arr in reg_rep.items()}

        # Create a function that takes a regular representation matrix and returns the corresponding
        # group element, using the reverse dictionary created above.
        def array_to_element(arr):
            return inv_reg_rep[array_to_tuple(arr)]

        return reg_rep, inv_reg_rep, element_to_array, array_to_element

    def verify_regular_representation(self, elem_to_arr, arr_to_elem):
        """Verifies that the regular representation satisfies the two requirements of it. This requires
        that the regular representation use dense matrices, NOT sparse matrices.
        """
        return ((self.identity == arr_to_elem(np.eye(self.order, dtype=int)))  # e == Vinv( identity_matrix )
                and
                # V(a) x V(b) == V(a * b)
                all([np.array_equal(np.dot(elem_to_arr(a), elem_to_arr(b)), elem_to_arr(self.op(a, b)))
                     for a in self
                     for b in self]))

    # ---------------------
    # Monoid Isomorphisms
    # ---------------------

    def make_element_mappings(self, other):
        """Returns a list of mappings (dictionaries) of this algebra's elements to all possible permutations
        of other's elements, where the identity of this algebra is always mapped to the identity of other."""
        if self.order == other.order:

            # Temporarily & non-destructively remove identities from lists of elements
            id0 = self.identity
            id1 = other.identity
            elems0copy = list(self.elements)
            elems1copy = list(other.elements)
            elems0copy.remove(id0)
            elems1copy.remove(id1)

            # Compute all possible mappings
            mappings = [dict(zip(elems0copy, perm)) for perm in it.permutations(elems1copy)]

            # Add the identities back in
            for mapping in mappings:
                mapping[id0] = id1
            return mappings
        else:
            raise ValueError(f"Algebras must be of the same order: {self.order} != {other.order}")


