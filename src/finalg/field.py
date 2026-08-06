# =========
#   Field
# =========

from finalg.group import Group
from finalg.ring import Ring
# from main import make_finite_algebra

def is_field(add_id, elements, table):
    """The elements of a Field, minus the additive identity, form a commutative Group
    under multiplication. This function takes the additive identity, the list of all
    elements, and a field's multiplication table as input, and returns the Group under
    multiplication, if it exists, otherwise it returns False.  If the proposed Field
    inputs are trivial (only one element and a 1x1 table) then False is returned.  That
    is, a trivial Field is not allowed."""
    from finalg.make_finite_algebra import make_finite_algebra
    if len(elements) == 1:
        return False
    else:
        mult = make_finite_algebra("tmp", "temporary", elements, table)
        # elems_copy = elements.copy()
        # elems_copy.remove(add_id)

        # elements is a tuple, which is immutable. But we want to remove the
        # additive identity element from it, so first turn elements into a list
        # then remove the additive identity, and finally turn the result back
        # into a tuple. Whew!
        elems_list = list(elements)
        elems_list.remove(add_id)
        elems_copy = tuple(elems_list)

        elems_copy_clo = mult.closure(elems_copy, True)  # Includes inverse elements
        if set(elems_copy) == set(elems_copy_clo):
            mult_sub = mult.subalgebra_from_elements(elems_copy)
            if isinstance(mult_sub, Group) and mult_sub.is_commutative():
                return mult_sub
            else:
                return False
        else:
            return False


class Field(Ring):
    """A Field is a Ring, where the elements, minus the additive identity, form a commutative Group
    under multiplication."""

    def __init__(self, name, description, elements, table, table2, check_inputs=True, mult_sub_grp=None,
                 conjugate_mapping=None):

        super().__init__(name, description, elements, table, table2, check_inputs, conjugate_mapping)

        # This is the abelian Group defined by the Ring elements, minus the additive identity,
        # under Ring multiplication
        self._mult_sub_grp = mult_sub_grp

        if check_inputs or mult_sub_grp is None:
            abelian_group = is_field(self.identity, self.elements, self.mult_table.table)
            if abelian_group:
                self._mult_sub_grp = abelian_group
            else:
                raise ValueError(f"Inputs do not support the construction of a Field.")

        self._mult_sub_grp.name = f"{self.name}_G"
        self._mult_sub_grp.description = f"Multiplicative abelian Group of {self.name}"

    def mult_abelian_subgroup(self):
        """Return the abelian Group defined by the Ring elements, minus the additive identity,
        under Ring multiplication."""
        return self._mult_sub_grp

    def mult_inv(self, element):
        """Return the multiplicative inverse of 'element', unless it's the additive identity
        element, in which case, return None."""
        if element == self.add_identity:
            return None
        else:
            return self._mult_sub_grp.inv(element)

    def div(self, x, y):
        """Return x/y, if y is not the additive identity; otherwise return None."""
        if y == self.add_identity:
            return None
        else:
            return self.mult(x, self.mult_inv(y))

    def element_to_power(self, elem, n, left_associative=True):
        """Overrides the Ring method by the same name, so that we use
        don't recreate the multiplicative Abelian subgroup already contained
        in the field.
        """
        mult_alg = self.mult_abelian_subgroup()
        return mult_alg.element_to_power(elem, n, left_associative)


