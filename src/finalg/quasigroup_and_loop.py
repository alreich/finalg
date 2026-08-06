# =============
#  Quasigroup
# =============

from finalg.magma import Magma


class Quasigroup(Magma):

    def __init__(self, name, description, elements, table):
        super().__init__(name, description, elements, table)


# =============
#     Loop
# =============

class Loop(Quasigroup):

    def __init__(self, name, description, elements, table):
        super().__init__(name, description, elements, table)


