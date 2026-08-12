from functools import reduce

# class FiniteOperator:
#     """A callable class that implements a binary operation based on a multiplication
#     table (i.e., Cayley table).  Although intended for use as the binary operation
#     of a finite algebra (e.g., Group operation), the implementation here can be called
#     with zero, one, two, or more arguments (similar to how arithmetic operators work in
#     Lisp).
#
#     If no arguments are provided, it will return the identity element if it exists;
#     otherwise it will return None.  e.g., op() ==> e | None
#
#     If only one argument is provided, it will check whether the argument is a valid
#     element of the algebra, and if so, return the same value, otherwise it will
#     raise an exception.  e.g., op(a) ==> a | ValueError
#
#     If two arguments are provided, it will return their 'product'.
#     e.g., op(a, b) ==> ab
#
#     If more than two arguments are provided, it will return their product by associating
#     left-to-right. e.g., op(a, b, c, d) = (((ab)c)d). The order of association is
#     only important for a Magma, because it is the only non-associative algebraic structure
#     supported here.
#     """
#
#     def __init__(self, elements, identity, table):
#         self._elements = elements
#         self._identity = identity
#         self._table = table
#
#     def __call__(self, *args):
#         return self._op(*args)
#
#     def _binary_operation(self, elem1, elem2):
#         """Returns the 'sum' of exactly two elements."""
#         row = self._elements.index(elem1)
#         col = self._elements.index(elem2)
#         index = self._table[row, col]
#         return self._elements[index]
#
#     def _op(self, *args):
#         if len(args) == 0:
#             return self._identity
#         elif len(args) == 1:
#             if args[0] in self._elements:
#                 return args[0]
#             else:
#                 raise ValueError(f"{args[0]} is not a valid element name")
#         elif len(args) == 2:
#             return self._binary_operation(args[0], args[1])
#         else:
#             return reduce(lambda a, b: self._op(a, b), args)

class FiniteOperator:
    """A callable class that implements a binary operation based on a multiplication
    table (i.e., Cayley table). Although intended for use as the binary operation
    of a finite algebra (e.g., Group operation), the implementation here can be called
    with zero, one, two, or more arguments (similar to how arithmetic operators work in
    Lisp).

    If no arguments are provided, it will return the identity element if it exists;
    otherwise it will return None.  e.g., op() ==> e | None

    If only one argument is provided, it will check whether the argument is a valid
    element of the algebra, and if so, return the same value, otherwise it will
    raise an exception.  e.g., op(a) ==> a | ValueError

    If two arguments are provided, it will return their 'product'.
    e.g., op(a, b) ==> ab

    If more than two arguments are provided, it will return their product by associating
    left-to-right. e.g., op(a, b, c, d) = (((ab)c)d). The order of association is
    only important for a Magma, because it is the only non-associative algebraic structure
    supported here.
    """

    def __init__(self, elements, identity, table):
        self._elements = elements
        self._identity = identity
        self._table = table
        # Name -> position lookup, built once, so _binary_operation doesn't have to
        # linearly re-scan self._elements (via list.index) on every single call. This
        # matters a lot for methods like Magma.closure(), which call the operation
        # O(k^2) times while growing a k-element subset -- each avoided O(n) scan adds
        # up fast once both k and n are non-trivial (e.g. subalgebra enumeration on A5).
        self._index = {e: i for i, e in enumerate(elements)}

    def __call__(self, *args):
        return self._op(*args)

    def _binary_operation(self, elem1, elem2):
        """Returns the 'sum' of exactly two elements."""
        try:
            row = self._index[elem1]
            col = self._index[elem2]
        except KeyError as exc:
            bad = elem1 if elem1 not in self._index else elem2
            raise ValueError(f"{bad} is not a valid element name") from exc
        index = self._table[row, col]
        return self._elements[index]

    def _op(self, *args):
        if len(args) == 0:
            return self._identity
        elif len(args) == 1:
            if args[0] in self._elements:
                return args[0]
            else:
                raise ValueError(f"{args[0]} is not a valid element name")
        elif len(args) == 2:
            return self._binary_operation(args[0], args[1])
        else:
            return reduce(lambda a, b: self._op(a, b), args)

