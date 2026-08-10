# ==============
# Element Class
# ==============

class Element:
    """Elements of the algebras here must be strings. This class turns a string element,
    along with its algebra, into an object with infix arithmetic operators, +, -, *, /, **,
    etc.
    """

    def __init__(self, name, algebra):

        self._algebra = algebra

        if isinstance(name, str):
            if name in self._algebra:
                self._name = name
            else:
                raise ValueError(f"name must be an element of algebra")
        else:
            raise ValueError(f"name must be a string")

    @property
    def name(self):
        """Return the name of this Element."""
        return self._name

    @property
    def algebra(self):
        """Return the algebra associated with this Element."""
        return self._algebra

    def __str__(self):
        return self._name

    def __repr__(self):
        return repr(self._name)

    def __add__(self, other):
        elem = self._algebra.op(self._name, other.name)
        return Element(elem, self._algebra)

    def __sub__(self, other):
        if hasattr(self._algebra, 'sub'):
            elem = self._algebra.sub(self._name, other.name)
            return Element(elem, self._algebra)
        else:
            raise ValueError(f"{self._algebra.name} does not support subtraction")

    def __neg__(self):
        if hasattr(self._algebra, 'inv'):
            elem = self._algebra.inv(self._name)
            return Element(elem, self._algebra)
        return None

    def __mul__(self, other):
        if hasattr(self._algebra, 'mult'):
            elem = self._algebra.mult(self._name, other.name)
            return Element(elem, self._algebra)
        else:
            raise ValueError(f"{self._algebra.name} does not support multiplication")

    def __truediv__(self, other):
        if hasattr(self._algebra, 'div'):
            elem = self._algebra.div(self._name, other.name)
            return Element(elem, self._algebra)
        else:
            raise ValueError(f"{self._algebra.name} does not support division")

    def __pow__(self, n):
        """See the documentation for the element_to_power method for either
        the Magma, Ring, or Field. This method calls one of those methods."""
        elem_to_pow_name = self._algebra.element_to_power(self._name, n)
        return Element(elem_to_pow_name, self._algebra)

    def __or__(self, other):
        """Conjugate the element with other element. That is, if a = self and b = other,
        then in Python notation, a ^ b returns inv(b) * a * b."""
        if hasattr(self._algebra, 'conjugate'):
            elem = self._algebra.conjugate(self._name, other.name)
            return Element(elem, self._algebra)
        else:
            raise ValueError(f"{self._algebra.name} does not support conjugation")

    def __key(self):
        return tuple([self._name, self._algebra.__hash__()])

    def __hash__(self):
        return hash(self.__key())

    def __eq__(self, other):
        if isinstance(other, Element):
            return self.__key() == other.__key() and self._algebra == other.algebra
        else:
            return NotImplemented


