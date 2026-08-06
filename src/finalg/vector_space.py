
from finalg.field import Field
from finalg.module import Module, module_sv_mult, module_dot_product, check_module_conditions


class VectorSpace(Module):
    """See https://abstract-algebra.readthedocs.io for the definition of a VectorSpace."""

    def __init__(self, name, description, field, group, operator):
        super().__init__(name, description, field, group, operator)
        if not isinstance(field, Field):
            raise ValueError(f"{field} must be a Field.")


class NDimensionalVectorSpace(VectorSpace):

    def __init__(self, field, n, check_input_conditions=True):
        name = f"{n}D-{field.name}"
        desc = f"{n}-dimensional Vector Space over {field.name}"
        self._dimensions = n

        # Group from the n-fold direct product of the Field with itself
        # group = field.power(n)
        group = field ** n

        super().__init__(name, desc, field, group, module_sv_mult(field))

        # Check input conditions, maybe
        if check_input_conditions:
            if not check_module_conditions(field, group, self.sv_mult):
                raise ValueError("Inputs don't meet required conditions.")

    @property
    def dimensions(self):
        """Returns the dimension of the VectorSpace's vectors."""
        return self._dimensions

    @property
    def origin(self):
        """Returns the origin element, a vector, of the VectorSpace."""
        return self.vector.identity

    def dot_product(self, u, v):
        """Computes and returns the dot-product of two VectorSpace vectors."""
        return module_dot_product(self, u, v)


