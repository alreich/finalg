# =================
#   FiniteAlgebra
# =================

from abc import ABC
import json

from finalg.cayley_table import CayleyTable, make_cayley_table
from finalg.element import Element

class FiniteAlgebra(ABC):
    """A top-level container class for functionality that is common to all finite algebras
    that only have one set of elements. THIS CLASS IS NOT INTENDED TO BE INSTANTIATED.
    """

    def __init__(self, name, description, elements, table):
        self.name = name
        self.description = description
        self._elements = tuple(elements)
        self._inverses = dict()

        # Set up the multiplication table
        if isinstance(table, CayleyTable):
            self._table = table
        else:
            self._table = make_cayley_table(table, elements)

        # If it exists, setup the algebra's identity element
        id_index = self._table.identity()
        if id_index is not None:
            self._identity = self.elements[id_index]
        else:
            self._identity = None

    def __contains__(self, element):
        return element in self._elements

    def __getitem__(self, index):
        return self._elements[index]

    def __repr__(self):
        nm = self.name
        desc = self.description
        elems = self._elements
        tbl = self._table.tolist()
        return f"{self.__class__.__name__}(\n'{nm}',\n'{desc}',\n{elems},\n{tbl}\n)"

    def __str__(self):
        return f"<{self.__class__.__name__}:{self.name}, ID:{id(self)}>"

    def __len__(self):
        """Same as the order of the algebra."""
        return len(self._elements)

    @property
    def elements(self):
        """Returns the algebra's element names (list of strings)."""
        return self._elements

    def element_map(self):
        """Instantiates an Element for each element name and returns a dictionary, where
         the element names are keys and corresponding Elements are the values. This method
         is used within the context manager, InfixNotation, to perform arithmetic using
         infix notation."""
        return {elem: Element(elem, self) for elem in self._elements}

    @property
    def table(self):
        """Returns the algebra's Cayley Table ('multiplication' table)."""
        return self._table

    @property
    def identity(self):
        """Returns the algebra's identity element if it exists; otherwise, it returns None."""
        return self._identity

    def has_identity(self):
        """A convenience function that returns True or False, depending on whether the algebra
        has an identity element."""
        if self._identity is None:
            return False
        else:
            return True

    @property
    def order(self):
        """Returns the order of the algebra."""
        return len(self._elements)

    def is_associative(self):
        """Returns True if the algebra is associative; returns False otherwise."""
        return self._table.is_associative()

    def is_commutative(self):
        """Returns True if the algebra is commutative; returns False otherwise."""
        return self._table.is_commutative()

    def is_abelian(self):
        """Returns True if the algebra is abelian; returns False otherwise."""
        return self.is_commutative()

    def has_cancellation(self, verbose=False):
        """Return True if, for every a & b in the algebra, there are unique x and y in the algebra
        such that ax=b and ya=b. Otherwise, return False. Set verbose to True to see intermediate
        calculations."""
        return self._table.has_cancellation(verbose)

    def has_inverses(self):
        """Returns True if every element in the algebra has an inverse that is also in the algebra;
        returns False otherwise."""
        return self._table.has_inverses()

    def create_inverse_lookup_dict(self):
        """Returns a dictionary that maps each of the algebra's elements to its inverse element."""
        if self.has_identity():
            return self._table.inverse_lookup_dict(self._table.identity(), self._elements)
        else:
            return None

    def inv(self, element):
        """Return the inverse of an element"""
        if self.has_inverses():
            if element in self._inverses:
                return self._inverses[element]
            else:
                return None
        else:
            return None

    def to_dict(self, include_classname=False):
        """Returns a dictionary that represents the algebra.  The dictionary
        can be fed back into make_finite_algebra, and it will return a copy of
        this algebra."""
        result = {'name': self.name,
                  'description': self.description,
                  'elements': self._elements,
                  'table': self._table.tolist()
                  }
        # If self is a Ring:
        if hasattr(self, 'mult_table'):
            result['table2'] = self.mult_table.tolist()
            if hasattr(self, 'conjugates'):
                if self.conjugates() is not None:
                    result['conj_map'] = self.conjugates()
        if include_classname:
            result['type'] = self.__class__.__name__
        return result

    def dumps(self):
        """Returns a JSON string that represents the algebra."""
        return json.dumps(self.to_dict())

    def dump(self, json_filename):
        """Writes the algebra to the given filename in JSON format."""
        with open(json_filename, 'w') as fout:
            json.dump(self.to_dict(), fout)
