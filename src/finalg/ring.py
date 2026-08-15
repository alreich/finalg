# ========
#   Ring
# ========

import numpy as np
import itertools as it
from re import finditer
from sympy.ntheory import isprime
from pprint import pprint

from finalg.utilities import delete_row_col, yes_or_no
from finalg.finite_operator import FiniteOperator
# from main import make_finite_algebra
from finalg.cayley_table import CayleyTable, make_cayley_table
from finalg.group import Group


class Ring(Group):
    """A Ring is a commutative Group with an 'addition' operator, along with an
    associative 'multiplication' operator, where multiplication distributes over
    addition.  The operator inherited from Group becomes 'addition', while
    'multiplication' is defined by a second Cayley table, table2."""

    def __init__(self, name, description, elements, table, table2, check_inputs=True,
                 conjugate_mapping=None):

        super().__init__(name, description, elements, table, check_inputs)

        self._conjugates = conjugate_mapping

        if isinstance(table2, CayleyTable):
            self._ring_mult_table = table2
        else:
            self._ring_mult_table = make_cayley_table(table2, elements)

        # If it exists, set up the Ring's multiplicative identity element
        mult_id_index = self._ring_mult_table.identity()
        if mult_id_index is not None:
            self._mult_identity = self.elements[mult_id_index]
        else:
            self._mult_identity = None

        self._ring_mult = FiniteOperator(self.elements, self._mult_identity, self._ring_mult_table)

        if check_inputs:
            if not super().is_commutative():
                raise ValueError(f"CHECK INPUTS: The ring addition operation is not commutative. {self}")
            if not self._ring_mult_table.is_associative():
                raise ValueError(f"CHECK INPUTS: The ring multiplication operation is not associative. {self}")
            if not self._ring_mult_table.distributes_over(self.table):
                raise ValueError(f"CHECK INPUTS: Multiplication does not distribute over addition. {self}")

    def __repr__(self):
        # nm, desc, elems, tbl, tbl2, conjmap = unpack_components(self)
        nm = self.name
        desc = self.description
        elems = self.elements
        tbl = self.table.tolist()
        tbl2 = self.mult_table.tolist()
        conjmap = self.conjugates()
        if conjmap is not None:
            return f"{self.__class__.__name__}(\n'{nm}',\n'{desc}',\n{elems},\n{tbl},\n{tbl2},\n{conjmap}\n)"
        else:
            return f"{self.__class__.__name__}(\n'{nm}',\n'{desc}',\n{elems},\n{tbl},\n{tbl2}\n)"

    def __mul__(self, other):  # Direct Product of two Rings
        """Return direct product of this Ring with the `other` Ring."""
        from finalg.make_finite_algebra import make_finite_algebra
        if not isinstance(other, Ring):
            raise ValueError(f"{other.name} must be a Ring")
        dp_name = f"{self.name}_x_{other.name}"
        dp_description = "Direct product of " + self.name + " & " + other.name
        dp_element_names = list(it.product(self.elements, other.elements))  # Cross product
        dp_add_table = list()
        dp_mul_table = list()
        for a in dp_element_names:
            dp_add_table_row = list()  # Start new rows in the add and mult tables
            dp_mul_table_row = list()
            for b in dp_element_names:
                dp_add_table_row.append(dp_element_names.index((self.add(a[0], b[0]),
                                                                other.add(a[1], b[1]))))
                dp_mul_table_row.append(dp_element_names.index((self.mult(a[0], b[0]),
                                                                other.mult(a[1], b[1]))))
            dp_add_table.append(dp_add_table_row)  # Append the new rows to each table
            dp_mul_table.append(dp_mul_table_row)
        return make_finite_algebra(dp_name,
                                   dp_description,
                                   list([f"{elem[0]}{self.direct_product_delimiter()}{elem[1]}"
                                         for elem in dp_element_names]),
                                   dp_add_table,
                                   dp_mul_table)

    def _key(self):
        return tuple([self.elements,
                      tuple(tuple(row) for row in self.table.tolist()),
                      tuple(tuple(row) for row in self._ring_mult_table.tolist())])

    def __hash__(self):
        return hash(self._key())

    def __eq__(self, other):
        if isinstance(other, Ring):
            return self._key() == other._key()
        else:
            return NotImplemented

    @property
    def add_identity(self):
        """Returns the additive identity element"""
        return self.identity

    @property
    def zero(self):
        """Another way to get the additive identity element"""
        return self.identity

    @property
    def mult_identity(self):
        """Returns the multiplicative identity element, if it exists.
        If it doesn't exist, then None is returned."""
        return self._mult_identity

    @property
    def one(self):
        """Another way to get the multiplicative identity element"""
        return self._mult_identity

    @property
    def minus_one(self):
        """Return the Ring's 'minus one' element. That is, the additive
        inverse of its multiplicative identity element."""
        return self.inv(self.mult_identity)

    def has_mult_identity(self):
        """A convenience function that returns True or False, depending on whether the algebra
        has a multiplicative identity element, in addition to its additive identity element."""
        if self.mult_identity is not None:
            return True
        else:
            return False

    @property
    def add_table(self):
        """Returns the CayleyTable for addition."""
        return self.table

    @property
    def mult_table(self):
        """Returns the CayleyTable for multiplication"""
        return self._ring_mult_table

    def add(self, *args):
        """Use the inherited group operator as the ring's addition operator."""
        return self.op(*args)

    def mult(self, *args):
        """Ring multiplication, based on the second table."""
        return self._ring_mult(*args)

    def element_to_power(self, elem, n, left_associative=True):
        """Overrides the Magma method by the same name, so that we use
        the multiplication operation of the Ring to raise an element to
        a power.
        """
        mult_alg = self.extract_multiplicative_algebra()
        return mult_alg.element_to_power(elem, n, left_associative)

    def mult_is_commutative(self):
        """By definition, Ring addition is commutative, but Ring multiplication only needs to be
        associative.  This method tells us whether multiplication is commutative for this Ring."""
        return self.mult_table.is_commutative()

    def extract_additive_algebra(self):
        """A Ring's elements over addition, alone, should be a commutative Group.  This function
        returns that Group."""
        from finalg.make_finite_algebra import make_finite_algebra
        nm = f"{self.name}.Add"
        desc = f"Additive-only portion of {self.name}"
        return make_finite_algebra(nm, desc, self.elements, self.table.table)

    def extract_multiplicative_algebra(self):
        """A Ring's elements over multiplication, alone, should be a Semigroup.  This function
        returns that Semigroup."""
        from finalg.make_finite_algebra import make_finite_algebra
        nm = f"{self.name}.Mult"
        desc = f"Multiplicative-only portion of {self.name}"
        return make_finite_algebra(nm, desc, self.elements, self.mult_table.table)

    # TODO: Write a method that returns non-zero pairs whose product is zero
    def zero_divisors(self):
        """Return the Ring's zero divisors. i.e., if neither a nor b are 0, but a*b == 0, then
        a and b are zero divisors."""

        # Get the index of the additive identity element ("zero")
        zero_index = self.elements.index(self.zero)

        # Delete the zero element's row & column in the multiplication table.
        # (NOTE: This operation leaves the original mult. table unchanged.)
        mult_table_without_add_id = delete_row_col(self.mult_table.table, zero_index, zero_index)

        # Get the row & column indices where the product equals "zero" in the remaining table
        a, b = list(map(lambda x: set(x), np.where(mult_table_without_add_id == zero_index)))
        #
        # Return all elements corresponding to the union of the row & column indices
        return [self.elements[index + 1] for index in list(a | b)]

    def units(self, return_names=True, verbose=False):
        """Return a list of the Ring's units."""
        mult_alg = self.extract_multiplicative_algebra()
        if mult_alg.has_identity():
            return mult_alg.units(return_names)
        else:
            if verbose:
                print(f"There is no multiplicative identity element.")
            return None

    def commutator(self, a, b):
        """Return [a, b] = (a * b) - (b * a), the ring commutator of a & b"""
        return self.sub(self.mult(a, b), self.mult(b, a))

    def element_pairs_where_product_equals(self, elem_name):
        """Return all pairs of elements where the product is equal to elem_name.
        """
        return self._element_pairs_where_table_equals(self.mult_table, elem_name)

    def zero_divisor_pairs(self):
        """Return a list of ordered pairs of elements, neither one of which is zero,
        but whose product is zero.
        """
        zero_product_pairs = self.element_pairs_where_product_equals(self.identity)
        # return [pair for pair in zero_product_pairs if not self.identity in pair]
        return [pair for pair in zero_product_pairs if self.identity not in pair]

    def square_root_mapping(self):
        """Return a dictionary where the keys are ring's elements and the values
        are the ring elements' square roots. Some elements may have no square
        roots, and some may have one or more square roots."""
        # The indices of elements with square roots are on the
        # diagonal of the multiplicative Cayley table.
        diag = self.mult_table.table.diagonal().tolist()
        elems_with_sqr_roots = set([self[elem] for elem in diag])
        # Create a dict with the necessary keys and empty lists
        result = {key: list() for key in elems_with_sqr_roots}
        for index in range(self.order):  # Populate the dict's empty lists
            key = self[diag[index]]
            val = self[index]
            result[key].append(val)
        return result

    def square_roots(self, elem_name):
        """Return a list of the square roots of elem_name. If the list is empty, there are none."""
        mapping = self.square_root_mapping()
        sqr_roots = list()
        if elem_name in mapping:
            sqr_roots = mapping[elem_name]
        return sqr_roots

    def about(self, max_size=12, max_gens=2, use_table_names=False, show_tables=True, show_elements=True,
              show_conjugates=False, show_generators=False):
        """Print information about the Ring."""
        super().about(max_size, max_gens, use_table_names, show_tables, show_elements, show_generators)

        if self.mult_identity is not None:
            print(f"Mult. Identity: {repr(self.mult_identity)}")
        else:
            print(f"Mult. Identity: None")

        print(f"Mult. Commutative? {yes_or_no(self.mult_is_commutative())}")

        zero_divisors = self.zero_divisors()
        if len(zero_divisors) == 0:
            print("Zero Divisors: None")
        else:
            print(f"Zero Divisors: {zero_divisors}")

        size = len(self.elements)

        if show_tables:
            if size <= max_size:
                if use_table_names:
                    print(f"Multiplicative Cayley Table (showing names):")
                    pprint(self._ring_mult_table.to_list_with_names(self.elements))
                else:
                    print(f"Multiplicative Cayley Table (showing indices):")
                    pprint(self.mult_table.tolist())
            else:
                print(f"{self.__class__.__name__} order is {size} > {max_size}, so the mult. table is not printed.")

        # The conjugate mapping is only printed out if it exists, and we want to see it.
        if show_conjugates and self.conjugates() is not None:
            print(f"Conjugate Mapping: {self.conjugates()}")

        return None

    # ========================================================================
    # The following methods implement the Cayley-Dickson construction/algebra
    # ========================================================================

    def sqr(self):  # My original version of the Cayley-Dickson construction/algebra
        """The Cayley-Dickson construction/algebra. Returns the direct product of this Ring
        with itself, where multiplication is defined as: (a, b) * (c, d) = (ac - bd, ad + bc)
        """
        from finalg.make_finite_algebra import make_finite_algebra
        dp_name = f"{self.name}_SQR"
        dp_description = "Direct product of " + self.name + " with itself using complex multiplication"
        dp_element_names = list(it.product(self.elements, self.elements))  # Cross product
        dp_add_table = list()
        dp_mul_table = list()
        for a in dp_element_names:
            dp_add_table_row = list()  # Start new rows in the add and mult tables
            dp_mul_table_row = list()
            for b in dp_element_names:
                dp_add_table_row.append(dp_element_names.index((self.add(a[0], b[0]),
                                                                self.add(a[1], b[1]))))
                # (a[0], a[1]) * (b[0], b[1])
                #     = (a[0]b[0] - a[1]b[1], a[0]b[1] + a[1]b[0])
                dp_mul_table_row.append(dp_element_names.index(((self.sub(self.mult(a[0], b[0]),
                                                                          self.mult(a[1], b[1]))),
                                                                (self.add(self.mult(a[0], b[1]),
                                                                          self.mult(a[1], b[0]))))))
            dp_add_table.append(dp_add_table_row)  # Append the new rows to each table
            dp_mul_table.append(dp_mul_table_row)
        return make_finite_algebra(dp_name,
                                   dp_description,
                                   list([f"{elem[0]}{self.direct_product_delimiter()}{elem[1]}"
                                         for elem in dp_element_names]),
                                   dp_add_table,
                                   dp_mul_table)

    def split_element(self, element):
        """If the element is a compound element created by a direct product or the Cayley-Dickson
        construction, then it contains at least one delimiter (e.g., '1:2' or '1:2:3:4').
        This method splits the element at the middle delimiter and returns the two pieces
        (e.g., '1', '2' or '1:2', '3:4').  If the element is not a compound element, then
        it is returned unchanged.
        """
        delimiter = self.direct_product_delimiter()
        if delimiter in element:
            matches = list(finditer(delimiter, element))
            mid = matches[len(matches) // 2]
            return element[:mid.start()], element[mid.end():]
        else:
            return element

    def is_gaussian_prime(self, elem):
        """This method only works for elements of Rings or Fields created by the function,
        'generate_algebra_mod_n', or by a single application of the Ring method,
        'make_cayley_dickson_algebra', to the output of 'generate_algebra_mod_n'.
        That is, elements that look like 7 or 07:12.
        """
        delim = self.direct_product_delimiter()
        dimension = elem.count(delim)
        if dimension == 0:
            real = int(elem)
            imag = 0
        elif dimension == 1:
            a, b = elem.split(delim)
            real = int(a)
            imag = int(b)
        else:
            raise ValueError(f"The dimension of {elem} is too high.")
        if real == 0:
            return isprime(imag) and imag % 4 == 3
        elif imag == 0:
            return isprime(real) and real % 4 == 3
        return isprime(real ** 2 + imag ** 2)

    def scalar_mult(self, scalar_name, elem_name):
        """ Scalar multiplication. 'a' * 'c:d' = 'a*c:a*d'
        Example: scalar_mult('2', '1:2', F3) ==> '2:1'
        """
        delimiter = self.direct_product_delimiter()
        scalarx = delimiter.join([scalar_name, self.zero[0]])  # eg: '2' --> '2:0'
        return self.mult(scalarx, elem_name)

    def elem_conj(self, elem):
        """For use only when making a Cayley-Dickson algebra.
        """
        delim = self.direct_product_delimiter()
        if delim in elem:
            a, b = self.split_element(elem)
            return delim.join([self.conj(a), self.inv(b)])
        else:
            return elem

    def conjugates(self):
        """Return the dictionary that maps elements to their conjugate values.
        If it's None, then the element is its own conjugate."""
        return self._conjugates

    def conj(self, elem):
        """Given an element name, return the element name of its conjugate value."""
        if self._conjugates is None:
            return elem
        else:
            return self._conjugates[elem]

    def norm(self, elem):
        """Return the product of the input element and its conjugate."""
        return self.mult(elem, self.conj(elem))

    def make_cayley_dickson_algebra(self, mu=None, version=1):
        """Constructs the Cayley-Dickson algebra using this Ring or Field.

        Several different versions of multiplication are supported:
        version=1: (DEFAULT) No mu & no conjugation are used
        version=2: Definition in Schafer, 1966
        version=3: Definition in Schafer, 1954
        version=4: Definition in Baez, 2001.

        See the documentation on readthedocs for more information regarding versions.

        Versions 2 & 3 require a value for mu. If mu is None (the default), then mu
        will be automatically set to be the additive inverse of the Ring's
        multiplicative identity element (i.e., "-1") if it exists. If it does not
        exist, then an exception will be raised.
        """
        from finalg.make_finite_algebra import make_finite_algebra
        if mu is None:
            if self.has_mult_identity():
                mu = self.inv(self.one)  # The additive inverse of the multiplicative identity
            else:
                if version == 2 or version == 3:
                    raise ValueError(f"Without a mult. identity, version {version} requires a specific value for mu")
        else:  # mu is not None
            if mu in self.elements:
                if version == 1 or version == 4:
                    print(f"** Version {version} ignores the value of mu. (mu = {mu})")
                    mu = None
            else:
                raise ValueError(f"mu = {mu} is not an element of {self.name}.")

        if version == 1:
            vers = "no mu and no conjugation were used."
            name_suffix = "_CDA_2024"
        elif version == 2:
            vers = f"mu = {mu}, Schafer 1966 version."
            name_suffix = "_CDA_1966"
        elif version == 3:
            vers = f"mu = {mu}, Schafer 1954 version."
            name_suffix = "_CDA_1954"
        elif version == 4:
            vers = f"mu = {mu}, Baez 2001 version."
            name_suffix = "_CDA_2001"
        else:
            raise ValueError(f"{version} is not a valid version #. Use 1, 2, or 3.")

        name = f"{self.name}{name_suffix}"
        description = f"Cayley-Dickson algebra based on {self.name}, where {vers}"
        element_names = list(it.product(self.elements, self.elements))  # Cross product
        elems = [f"{elem[0]}{self.direct_product_delimiter()}{elem[1]}" for elem in element_names]

        # The conjugate mapping created here will be passed, at the end, to the CD algebra
        # output by this method.  On the other hand, the self.conj method calls in the
        # multiplication section, below, refer to the conjugates of this ring itself, not
        # the ring to be output.
        conj_elems = [self.elem_conj(elem) for elem in elems]
        conj_map = dict(zip(elems, conj_elems))

        add_table = list()
        mul_table = list()
        for x in element_names:
            a = x[0]
            b = x[1]

            # Start new rows in the addition and multiplication tables
            add_table_row = list()
            mul_table_row = list()
            for y in element_names:
                c = y[0]
                d = y[1]

                # ADDITION: (a, b) + (c, d) = (a + b, c + d)
                add_table_row.append(element_names.index((self.add(a, c),
                                                          self.add(b, d))))
                # MULTIPLICATION:
                if version == 1:
                    # No mu and no conjugation used:
                    # Multiplication: (a, b) x (c, d) = (a x c  -  b x d,  a x d  +  b x c)
                    mul_table_row.append(element_names.index(((self.sub(self.mult(a, c),
                                                                        self.mult(b, d))),
                                                              (self.add(self.mult(a, d),
                                                                        self.mult(b, c))))))
                elif version == 2:
                    # See [Schafer, 1966]
                    # Conjugation: a* = a and (u, v)* = (u*, -v) recursively
                    # Multiplication: (a, b) x (c, d) = (a x c  +  mu x d x b*,  a* x d  +  c x b)
                    mul_table_row.append(element_names.index(((self.add(self.mult(a, c),
                                                                        self.mult(mu, d, self.conj(b)))),
                                                              (self.add(self.mult(self.conj(a), d),
                                                                        self.mult(c, b))))))
                elif version == 3:
                    # See [Schafer, 1954]
                    # Multiplication: (a, b) x (c, d) = (a x c  +  mu x d* x b,  d x a  +  b x c*)
                    mul_table_row.append(element_names.index(((self.add(self.mult(a, c),
                                                                        self.mult(mu, self.conj(d), b))),
                                                              (self.add(self.mult(d, a),
                                                                        self.mult(b, self.conj(c)))))))
                elif version == 4:
                    # See [Baez 2001]
                    # Multiplication: (a, b) x (c, d) = (a x c  -  d x b*,  a* x d  +  c x b)
                    mul_table_row.append(element_names.index(((self.sub(self.mult(a, c),
                                                                        self.mult(d, self.conj(b)))),
                                                              (self.add(self.mult(self.conj(a), d),
                                                                        self.mult(c, b))))))
                else:
                    raise ValueError(f"What happened?!?! We should never see this message. Version == {version}")

            # Append the new rows to each table
            add_table.append(add_table_row)
            mul_table.append(mul_table_row)

        return make_finite_algebra(name,
                                   description,
                                   elems,
                                   add_table,
                                   mul_table,
                                   conj_map)


