# =============
#   Semigroup
# =============

from finalg.magma import Magma

class Semigroup(Magma):
    """A Semigroup is an associative Magma."""

    def __init__(self, name, description, elements, table, check_inputs=True):
        super().__init__(name, description, elements, table)
        if check_inputs:
            if not self.table.is_associative():
                raise ValueError("CHECK INPUTS: Table does not support associativity")

    def is_regular(self):
        """Returns True if for all elements, a, there exists an element, x, such that axa=a.
        Otherwise, False is returned."""
        return all([any([self.op(self.op(a, x), a) == a for x in self]) for a in self])

    def weak_inverses(self):
        """Returns a dictionary of weak inverses, where each key is one of the algebra's
        elements and its value is a list of its weak inverses.  If the algebra is
        regular, then there will be at least 1 weak inverse for each element. Otherwise,
        some elements may not have a weak inverse."""
        return {a: [x for x in self if self.op(self.op(a, x), a) == a] for a in self}


