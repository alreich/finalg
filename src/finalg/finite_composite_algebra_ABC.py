from abc import ABC

class FiniteCompositeAlgebra(ABC):
    """This class represents Finite Algebras that have more than one element list,
    such as VectorSpaces and Modules."""

    def __init__(self, name: str, description: str):
        self.name = name
        self.description = description


