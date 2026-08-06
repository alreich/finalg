
class InfixNotation:
    """Creates a context manager for doing finite algebra calculations using infix notation.

    EXAMPLE:

    s3 = generate_symmetric_group(3)
    with InfixNotation(s3) as f:
        print(f['(2, 1, 3)'] + f['(3, 2, 1)'])

    ==> (3, 1, 2)

    """

    def __init__(self, algebra):
        self.element_map = algebra.element_map()

    def __enter__(self):
        return self.element_map

    def __exit__(self, _type, value, traceback):
        pass


