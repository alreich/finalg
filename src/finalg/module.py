
from functools import reduce

from finalg.finite_composite_algebra_ABC import FiniteCompositeAlgebra
from finalg.group import Group
from finalg.ring import Ring
from finalg.utilities import yes_or_no


class Module(FiniteCompositeAlgebra):
    """See https://abstract-algebra.readthedocs.io for the definition of a Module"""

    def __init__(self, name, description, ring, group, operator):
        super().__init__(name, description)
        if not isinstance(ring, Ring):
            raise ValueError(f"{ring} is not a Ring.")
        if not isinstance(group, Group) and group.is_abelian():
            raise ValueError(f"{group} is not an abelian Group.")
        if not check_module_conditions(ring, group, operator):
            raise ValueError("Inputs don't meet requirements for a Module.")
        self.scalar = ring
        self.vector = group
        self.sv_mult = operator  # scalar-vector operator

    def __repr__(self):
        sname = self.scalar.name
        vname = self.vector.name
        cname = self.__class__.__name__
        return f"<{cname}:{self.name}, ID:{id(self)}, Scalars:{sname}, Vectors:{vname}>"

    def vector_add(self, v1, v2):
        """Return the sum of two vectors using the Group operation, op."""
        return self.vector.op(v1, v2)

    def about(self, max_size=12, max_gens=2, use_table_names=False, show_tables=True, show_elements=True,
              show_generators=False):
        """Print information about the Module or Vector Space."""
        print(f"\n{self.__class__.__name__}: {self.name}")
        print(f"Instance ID: {id(self)}")
        print(f"Description: {self.description}")
        print(f"\nSCALARS:")
        self.scalar.about(max_size, max_gens, use_table_names, show_tables, show_elements, show_generators)
        print(f"\nVECTORS:")
        self.vector.about(max_size, max_gens, use_table_names, show_tables, show_elements, show_generators)
        return None


class NDimensionalModule(Module):

    def __init__(self, ring, n, check_input_conditions=True):
        name = f"{n}D-{ring.name}"
        desc = f"{n}-dimensional Module over {ring.name}"
        self._dimensions = n

        # Group from the n-fold direct product of the Field with itself
        # group = ring.power(n)
        group = ring ** n

        super().__init__(name, desc, ring, group, module_sv_mult(ring))

        # Check input conditions, maybe
        if check_input_conditions:
            if not check_module_conditions(ring, group, self.sv_mult):
                raise ValueError("Inputs don't meet required conditions.")

    @property
    def dimensions(self):
        """Returns the dimension of the Module's vectors."""
        return self._dimensions

    @property
    def origin(self):
        """Returns the origin element, a vector, of the Module."""
        return self.vector.identity

    def dot_product(self, u, v):
        """Computes and returns the dot-product of two Module vectors."""
        return module_dot_product(self, u, v)


# TODO: move this to be inside NDimensionalModule
def module_sv_mult(ring):
    """Returns a function that scales a vector.  That is, a function that takes
    a scalar and a vector, and returns their product, also a vector."""
    delimiter = ring.direct_product_delimiter()

    # sv_mult(s, v) takes an element created from a direct product (e.g., v = "a:b:c"),
    # splits it into a list (e.g., ["a", "b", "c"]), then maps the multiplication of
    # another element, say "s", over the list (e.g., ["s" * "a", "s" * "b", "s" * "c"])
    # and then joins the list back together into a single string (e.g., "sa:sb:sc"),
    # where sa, sb, & sc represent the results of the multiplications.
    def sv_mult(s, v):
        """Scalar-Vector product function"""
        return delimiter.join([ring.mult(s, x) for x in v.split(delimiter)])

    return sv_mult


# TODO: move this to be inside NDimensionalModule
def module_dot_product(ring, vec1, vec2):
    """Returns a scalar (ring element) that represents the dot-product of the
    two input vectors."""
    delim = ring.scalar.direct_product_delimiter()
    return reduce(lambda a, b: ring.scalar.add(a, b),
                  map(lambda pair: ring.scalar.mult(*pair),
                      zip(vec1.split(delim), vec2.split(delim))))


def check_module_conditions(ring: Ring, group: Group, sv_mult, verbose=False):
    """Returns True if all four conditions required of a Module hold true,
    otherwise this function returns False."""

    check1 = check_scaling_by_one(ring, group, sv_mult, verbose)
    if verbose:
        print(f"* Scaling by 1 OK? {yes_or_no(check1)}")

    check2 = check_dist_of_scalars_over_vec_add(ring, group, sv_mult, verbose)
    if verbose:
        print(f"* Distributivity of scalars over vector addition OK? {yes_or_no(check2)}")

    check3 = check_dist_of_vec_over_scalar_add(ring, group, sv_mult, verbose)
    if verbose:
        print(f"* Distributivity of vectors over scalar addition OK? {yes_or_no(check3)}")

    check4 = check_associativity(ring, group, sv_mult, verbose)
    if verbose:
        print(f"* Scaling by 1 OK? {yes_or_no(check4)}")

    return check1 & check2 & check3 & check4


def check_scaling_by_one(ring, group, sv_mult, verbose=False):
    """Returns True if scaling by one holds true in all cases, otherwise False is Returned."""
    is_ok = True
    one = ring.one
    for v in group.elements:
        if v != sv_mult(one, v):
            is_ok = False
            if verbose:
                print(f"{one} x {v} = {sv_mult(one, v)}")
    return is_ok


def check_dist_of_scalars_over_vec_add(ring, group, sv_mult, verbose=False):
    """Returns True if distributivity of scalars over vector addition holds true in all cases,
    otherwise False is Returned."""
    is_ok = True
    for s in ring.elements:
        for v1 in group.elements:
            for v2 in group.elements:
                a = sv_mult(s, group.op(v1, v2))
                b = group.op(sv_mult(s, v1), sv_mult(s, v2))
                if a != b:
                    is_ok = False
                    if verbose:
                        print(f"{a} != {b}")
    return is_ok


def check_dist_of_vec_over_scalar_add(ring, group, sv_mult, verbose=False):
    """Returns True if distributivity of vectors over scalar addition holds true in all cases,
    otherwise False is Returned."""
    is_ok = True
    for s1 in ring.elements:
        for s2 in ring.elements:
            for v in group.elements:
                a = sv_mult(ring.add(s1, s2), v)
                b = group.op(sv_mult(s1, v), sv_mult(s2, v))
                if a != b:
                    is_ok = False
                    if verbose:
                        print(f"{a} != {b}")
    return is_ok


def check_associativity(ring, group, sv_mult, verbose=False):
    """Return True if the special associativity condition on scalars and vectors holds true,
    otherwise return False."""
    is_ok = True
    for s1 in ring.elements:
        for s2 in ring.elements:
            for v in group.elements:
                a = sv_mult(ring.add(s1, s2), v)
                b = group.op(sv_mult(s1, v), sv_mult(s2, v))
                if a != b:
                    is_ok = False
                    if verbose:
                        print(f"{a} != {b}")
    return is_ok




