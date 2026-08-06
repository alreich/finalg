Groups, Monoids, Semigroups, Quasigroups, Loops, & Magmas
=========================================================

This section provides numerous examples of finite algebra creation and
manipulation, specifically, for algebras with only one set of elements
and one binary operation: Groups, Monoids, Semigroups, Quasigroups,
Loops, and Magmas. See the previous section, “Definitions”, for
definitions of these algebraic structures.

Representation of a FiniteAlgebra with only one Operation
---------------------------------------------------------

Internally, a ``FiniteAlgebra`` can take several different forms. For
algebras that have only one set of elements and one binary operation,
such as Groups, Monoids, Semigroups, Quasigroups, Loops, & Magmas, the
internal representation is described below.

-  **name**: (``str``) A short name for the algebra;
-  **description**: (``str``) Any additional, useful information about
   the algebra;
-  **elements**: (``list`` of ``str``) Names of the algebras’s elements.
-  **table**: (``list`` of ``list`` of ``int``) The algebra’s
   multiplication table, where each list in the list represents a row of
   the table, and each integer represents the position of an element in
   ‘element_names’. For input and ouput, element names (``str``) may be
   used in the table, rather than integers, but integers are still used
   internally.

The following section describes the required table in more detail.

Cayley Table
------------

A binary operation over a set of order :math:`n`, can be represented by
a 2-dimensional, :math:`n \times n` table, called a *Cayley Table*. All
of the properties of a finite algebra can be derived from its Cayley
Table.

For example, the Cayley table for the `commutative, non-associative
Magma derived from the Rock (r), Paper (p), Scissors (s)
game <https://en.wikipedia.org/wiki/Commutative_magma>`__, is shown on
the left in the diagram below. Rock, Paper, Scissors is a type of
algebra known as *Magma*.

If the elements are assumed to be represented in an ordered list, such
as (r, p, s), then the Cayley table can be unambiquously represented,
without row or column headers, as a list of integers that represent the
zero-based position of the elements in the ordered list, as shown on the
right in the diagram below.

.. figure:: ../docs/_static/Cayley_Table_Diagram.png
   :alt: The Cayley Table Representation

   The Cayley Table Representation

Importing the finalg module
---------------------------

Currently, when the ``finalg`` module is imported, a list of example
algebras is automatically retrieved and stored in the variable,
``examples``, and a table of contents (TOC) is printed out. Any specific
algebra in the list can be retrieved using its index, as shown in the
TOC.

.. code:: ipython3

    >>> import finalg as fa


.. parsed-literal::

    
    Loading examples...
    
    To retrieve an example, use finalg.examples[INDEX]
    
    ======================================================================
                               Example Algebras
    ----------------------------------------------------------------------
      19 example algebras are available.
      The INDEX is the first number on each line below:
    ----------------------------------------------------------------------
    0: A4 -- Alternating group on 4 letters (AKA Tetrahedral group)
    1: D3 -- https://en.wikipedia.org/wiki/Dihedral_group_of_order_6
    2: D4 -- Dihedral group on four vertices
    3: Pinter29 -- Non-abelian group, p.29, 'A Book of Abstract Algebra' by Charles C. Pinter
    4: RPS -- Rock, Paper, Scissors Magma
    5: S3 -- Symmetric group on 3 letters
    6: S3X -- Another version of the symmetric group on 3 letters
    7: V4 -- Klein-4 group
    8: Z4 -- Cyclic group of order 4
    9: F4 -- Field with 4 elements (from Wikipedia)
    10: mag_id -- Magma with Identity
    11: Example 1.4.1 -- See: Groupoids and Smarandache Groupoids by W. B. Vasantha Kandasamy
    12: Ex6 -- Example 6: http://www-groups.mcs.st-andrews.ac.uk/~john/MT4517/Lectures/L3.html
    13: Q8 -- Quaternion Group
    14: SD16 -- Semidihedral group of order 16
    15: A5 -- Alternating group on 5 letters
    16: F2 -- Field with 2 elements from paper: 236w06fields.pdf
    17: Latin_Sqr -- Latin Square. A division algebra (AKA Quasigroup)
    18: IP_Loop -- IP loop of small order
    ======================================================================


The Rock-Paper-Scissors algebra is one of the example algebras:

.. code:: ipython3

    >>> rps = fa.examples[4]
    >>> rps




.. parsed-literal::

    Magma(
    'RPS',
    'Rock, Paper, Scissors Magma',
    ('r', 'p', 's'),
    [[0, 1, 0], [1, 1, 2], [0, 2, 2]]
    )



All algebras have an ``about`` method, that will provide information
about the algebra’s structure and properties:

.. code:: ipython3

    >>> rps.about()


.. parsed-literal::

    
    ** Magma **
    Name: RPS
    Instance ID: 5215592736
    Description: Rock, Paper, Scissors Magma
    Order: 3
    Identity: None
    Associative? No
    Commutative? Yes
    Cyclic?: No
    Elements: ('r', 'p', 's')
    Has Cancellation? No
    Has Inverses? No
    Cayley Table (showing indices):
    [[0, 1, 0], [1, 1, 2], [0, 2, 2]]


A more detailed look at Rock-Paper-Scissors is coming up, below, along
with more information on Magmas.

The following sections describe the preferred algebra construction
process.

Algebra Constuction Examples
----------------------------

In a nutshell, use the function, ``make_finite_algebra`` for all algebra
construction.

Although individual algebras (Magma, Semigroup, etc.) have their own
individual constructors, the **recommended** way to construct an algebra
is to use the function, ``make_finite_algebra``, using one of the
following three approaches to inputs:

1. Enter **individual values** corresponding to the quantities in its
   Internal Representation, described above.
2. Enter a **Python dictionary** (``dict``), with keys and values
   corresponding to the individual values, described above.
3. Enter the **path to a JSON file** (``str``) that corresponds to the
   dictionary, described above. The examples use this method.

``make_finite_algebra`` uses the input table to determine what type of
algebra it supports and returns the appropriate algebra.

In the following examples, the only algebra constructor used is
``make_finite_algebra``.

Group
~~~~~

We’ll start the examples in the middle of the hierarchy of algebras,
with Groups.

The element names in a finite algebra’s ordered element list (actually
an immutable tuple) are **always** represented as strings; and, although
a Cayley table can be entered (and displayed) using the same string
names, within a Cayley table they are represented and displayed using a
2-dimensional, square array of integers that denote the positions of
element names in the element list. To see this, look at the following
input table vs. output table.

.. code:: ipython3

    >>> z3 = fa.make_finite_algebra('Z3',
                                    'Cyclic group of order 3',
                                    ['e', 'a', 'a^2'],  # A list of elements can be entered, but will be stored as a tuple
                                    [[ 'e' ,  'a' , 'a^2'],
                                     [ 'a' , 'a^2',  'e' ],
                                     ['a^2',  'e' ,  'a' ]]
                                   )
    >>> z3




.. parsed-literal::

    Group(
    'Z3',
    'Cyclic group of order 3',
    ('e', 'a', 'a^2'),
    [[0, 1, 2], [1, 2, 0], [2, 0, 1]]
    )



Following Python convention, the representation (``__repr__``) of an
algebra can be used to construct the same algebra, via copy-and-paste.

On the other hand, printing an algebra converts it to a compact string
representation that contains its class name, algebra name, and the
unique ID of the algebra instance:

.. code:: ipython3

    >>> print(z3)


.. parsed-literal::

    <Group:Z3, ID:5182807344>


The ``about`` method prints information about an algebra. Set
``use_element_names`` to ``True`` to see the Cayley table printed using
element names (``str``) rather than element positions (``int``).

.. code:: ipython3

    >>> z3.about(use_table_names=True)


.. parsed-literal::

    
    ** Group **
    Name: Z3
    Instance ID: 5182807344
    Description: Cyclic group of order 3
    Order: 3
    Identity: 'e'
    Commutative? Yes
    Cyclic?: Yes
    Generators: ['a', 'a^2']
    Elements:
       Index   Name   Inverse  Order
          0     'e'     'e'       1
          1     'a'   'a^2'       3
          2   'a^2'     'a'       3
    Cayley Table (showing names):
    [['e', 'a', 'a^2'], ['a', 'a^2', 'e'], ['a^2', 'e', 'a']]




.. parsed-literal::

    '<Group:Z3, ID:5182807344>'



Group Properties
~~~~~~~~~~~~~~~~

.. code:: ipython3

    >>> z3.is_associative()  # Only Magmas are non-associative




.. parsed-literal::

    True



.. code:: ipython3

    >>> z3.is_commutative()




.. parsed-literal::

    True



.. code:: ipython3

    >>> z3.is_abelian()




.. parsed-literal::

    True



The ``identity`` method (property) returns the algebra’s identity
element, if it exists. This method is implemented as a property, so no
trailing parentheses (“()”) are required. See the API documentation to
see what other methods are implemented as properties.

If the identity doesn’t exist, then ``None`` is returned.

.. code:: ipython3

    >>> z3.identity




.. parsed-literal::

    'e'



.. code:: ipython3

    >>> z3.inv('a')  # Get an element's inverse, if it exists




.. parsed-literal::

    'a^2'



Internal to algebras, tables are stored as instances of the
``CayleyTable`` class. Under normal usage, there should be no need to
work directly with a ``CayleyTable``.

.. code:: ipython3

    >>> z3.table




.. parsed-literal::

    CayleyTable([[0, 1, 2], [1, 2, 0], [2, 0, 1]])



Binary Operation
~~~~~~~~~~~~~~~~

The binary operation, implicitely defined by a Cayley table, is made
explicit by an algebra’s ``op`` method. Obviously, ``op`` should be able
to take two algebraic elements as input, and it does, but it can also
take any number of inputs from none to as many as needed. This was done
as a convenience, because the binary operation here is implemented as a
function, e.g., ``op('a', 'a^2')``, rather than as an infix operator.
This way, expressions such as :math:`b \cdot a \cdot b^{-1}`, can be
computed by ``op('b', 'a', g.inv('b'))``, instead of
``op(op('b', 'a'), g.inv('b'))``.

Note, however, that for more than two arguments, the binary operation is
performed left-to-right. That is,
:math:`b \cdot a \cdot b^{-1} \equiv (b \cdot a) \cdot b^{-1}`. This
matters for non-associative algebras, such as Magmas, Quasigroups, and
Loops.

For the example Group, Z3, created above, we have that
:math:`a \circ a = a^2` and this is verified below.

.. code:: ipython3

    >>> z3.op('a', 'a')




.. parsed-literal::

    'a^2'



Also, the following holds,
:math:`a \circ a \circ a = a \circ a^2 = a^2 \circ a = e`, as shown
below

.. code:: ipython3

    >>> z3.op('a', 'a', 'a') == z3.op('a', 'a^2') == z3.op('a^2', 'a') == 'e'




.. parsed-literal::

    True



The operations depicted in the two examples, above, can also be
accomplished using the context manager, **InfixNotation**, using with-as
syntax, as shown below. The target, **z**, of the context manager is a
dictionary where the keys are the names (str) of the algebra’s elements,
and the values are **Element** objects that can be used in infix-base
arithmetic expressions. So, for example, ``z['a']`` retrieves the
Element object that corresponds to the element named ‘a’.

.. code:: ipython3

    >>> with fa.InfixNotation(z3) as z:
    >>>     print(z['a'] + z['a'])
    >>>     print(z['a'] + z['a^2'])


.. parsed-literal::

    a^2
    e


If only one argument is given to an algebra’s binary operation, ``op``,
then that argument is simply returned; unless it is not a valid element
of the algebra, in which case an exception is raised.

.. code:: ipython3

    >>> z3.op('a')




.. parsed-literal::

    'a'



With zero arguments, ``op`` returns the identity, if it exists.

.. code:: ipython3

    >>> z3.op()




.. parsed-literal::

    'e'



Note, however, that ``op`` can only be used with elements (``str``) that
are members of an algebra’s element list. And, since ‘a^3’ is not a
string in Z3’s element list, it cannot be used in ``op``.

.. code:: ipython3

    >>> try:
    >>>     z3.op('a^3')
    >>> except Exception as exc:
    >>>     print(exc)


.. parsed-literal::

    a^3 is not a valid element name


However, if ‘a’ raised the power 3 is desired, then the following works:

.. code:: ipython3

    >>> with fa.InfixNotation(z3) as z:
    >>>    print(z['a'] ** 3)


.. parsed-literal::

    e


As does the following,

.. code:: ipython3

    >>> z3.element_to_power('a', 3)




.. parsed-literal::

    'e'



“Subtraction” in Groups
~~~~~~~~~~~~~~~~~~~~~~~

The method, ``sub``, is a convenience method for computing
“:math:`x - y`”, that is, :math:`x \circ y^{-1}` where
:math:`x, y \in S` and :math:`\langle S, \circ \rangle` is a Group.

.. code:: ipython3

    >>> x = 'a'
    >>> y = 'a^2'
    >>> print(f"For example, \"{x} - {y}\" = {x} * inv({y}) = {x} * {z3.inv(y)} = {z3.op(x, z3.inv(y))}")


.. parsed-literal::

    For example, "a - a^2" = a * inv(a^2) = a * a = a^2


Or, more succinctly using ``sub``:

.. code:: ipython3

    >>> z3.sub(x, y)




.. parsed-literal::

    'a^2'



Or, using infix notation via the Algebra context manager:

.. code:: ipython3

    >>> with fa.InfixNotation(z3) as z:
    >>>     print(z['a'] - z['a^2'])


.. parsed-literal::

    a^2


Magma
~~~~~

A Magma is the most fundamental type of algebra that we can instantiate
with ``make_finite_algebra``.

**Magma** – :math:`\langle S, \circ \rangle`, where :math:`S` is a set
and :math:`\circ` is a binary operation, :math:`\circ: S \times S \to S`

**Example: Rock-Paper-Scissors**

-  paper covers rock
-  rock crushes scissors
-  scissors cuts paper

Expressing this in algebraic form (see
https://en.wikipedia.org/wiki/Commutative_magma), where p *beats* r, and
r *beats* s, and s *beats* p, we have:

-  :math:`\langle S, \circ \rangle`, where :math:`S = \{r,p,s\}`
-  For all :math:`x, y \in S`, if :math:`x` *beats* :math:`y`, then
   :math:`x \circ y = y \circ x = x`
-  Also, for all :math:`x \in S`, :math:`x \circ x = x`

From the rule in the second bullet, above, this algebra is obviously
commutative.

.. code:: ipython3

    >>> rps = fa.make_finite_algebra('RPS',
                                     'Rock, Paper, Scissors Magma',
                                     ['r', 'p', 's'],
                                     [['r', 'p', 'r'],
                                      ['p', 'p', 's'],
                                      ['r', 's', 's']])
    
    >>> rps.about()


.. parsed-literal::

    
    ** Magma **
    Name: RPS
    Instance ID: 5215433296
    Description: Rock, Paper, Scissors Magma
    Order: 3
    Identity: None
    Associative? No
    Commutative? Yes
    Cyclic?: No
    Elements: ('r', 'p', 's')
    Has Cancellation? No
    Has Inverses? No
    Cayley Table (showing indices):
    [[0, 1, 0], [1, 1, 2], [0, 2, 2]]


Paper beats Rock:

.. code:: ipython3

    >>> rps.op('r', 'p')




.. parsed-literal::

    'p'



.. code:: ipython3

    >>> if rps.identity is None:
    >>>     print("RPS does not have an identity element")


.. parsed-literal::

    RPS does not have an identity element


For convenience, the method, ``has_identity``, returns True or False,
depending on whether an algebra has an identity.

.. code:: ipython3

    >>> rps.has_identity()




.. parsed-literal::

    False



The next section demonstrates that a Magma can have an identity element,
as long as the Magma is not associative, otherwise
``make_finite_algebra`` would output a Monoid.

A convention often used with abstract algebras is to denote an identity
element with the letter “e”. That is done in the example below, but the
function, ``make_finite_algebra``, ignores that convention and derives
the identity element, if it exists, from the table properties alone.

Magma with Identity Element
~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code:: ipython3

    >>> mag = fa.make_finite_algebra('Whatever',
                                     'Magma with Identity',
                                     ['e', 'a', 'b'],
                                     [['e', 'a', 'b'],
                                      ['a', 'e', 'a'],
                                      ['b', 'b', 'a']])
    
    >>> mag.about()


.. parsed-literal::

    
    ** Magma **
    Name: Whatever
    Instance ID: 5215720272
    Description: Magma with Identity
    Order: 3
    Identity: e
    Associative? No
    Commutative? No
    Cyclic?: Yes
    Generators: ['b']
    Elements: ('e', 'a', 'b')
    Has Cancellation? No
    Has Inverses? No
    Cayley Table (showing indices):
    [[0, 1, 2], [1, 0, 1], [2, 2, 1]]


Quasigroup
~~~~~~~~~~

Magmas can sometimes possess the property of **cancellation** that
appears similar to the ability to perform *division*, but is really far
from it. Nevertheless, such algebras are said to **have division**.

A Magma, :math:`M = \langle S, \circ \rangle` has the **Cancellation
Property** if :math:`\forall a,b \in S, \exists !\ x,y \in S` such that
:math:`a \circ x = b` and :math:`y \circ a = b`.

(Note that, in the definition, above, :math:`x` and :math:`y` must exist
*uniquely*.)

A Magma with the cancellation property is called a **Quasigroup** or
**division Magma**.

**Quasigroup Example**

The following description and example of a Quasigroup are from
Wikipedia: *“The multiplication table of a finite quasigroup is a Latin
square: an n × n table filled with n different symbols in such a way
that each symbol occurs exactly once in each row and exactly once in
each column.”* (https://en.wikipedia.org/wiki/Quasigroup#Latin_squares)

.. code:: ipython3

    >>> latin_sqr = fa.make_finite_algebra('Latin_Sqr',
                                           'Latin Square. A division algebra (AKA Quasigroup)',
                                           ('a0', 'a1', 'a2', 'a3', 'a4', 'a5', 'a6', 'a7', 'a8', 'a9'),
                                           [[0, 4, 8, 2, 3, 9, 6, 7, 1, 5], 
                                            [3, 6, 2, 8, 7, 1, 9, 5, 0, 4], 
                                            [8, 9, 3, 1, 0, 6, 4, 2, 5, 7], 
                                            [1, 7, 6, 5, 4, 8, 0, 3, 2, 9], 
                                            [2, 1, 9, 0, 6, 7, 5, 8, 4, 3], 
                                            [5, 2, 7, 4, 9, 3, 1, 0, 8, 6], 
                                            [4, 3, 0, 6, 1, 5, 2, 9, 7, 8], 
                                            [9, 8, 5, 7, 2, 0, 3, 4, 6, 1], 
                                            [7, 0, 1, 9, 5, 4, 8, 6, 3, 2], 
                                            [6, 5, 4, 3, 8, 2, 7, 1, 9, 0]])
    
    >>> latin_sqr.about()


.. parsed-literal::

    
    ** Quasigroup **
    Name: Latin_Sqr
    Instance ID: 5215433616
    Description: Latin Square. A division algebra (AKA Quasigroup)
    Order: 10
    Identity: None
    Associative? No
    Commutative? No
    Cyclic?: Yes
    Generators: ['a9', 'a5', 'a7', 'a8', 'a4', 'a3', 'a6', 'a1', 'a2']
    Elements: ('a0', 'a1', 'a2', 'a3', 'a4', 'a5', 'a6', 'a7', 'a8', 'a9')
    Has Cancellation? Yes
    Has Inverses? No
    Cayley Table (showing indices):
    [[0, 4, 8, 2, 3, 9, 6, 7, 1, 5],
     [3, 6, 2, 8, 7, 1, 9, 5, 0, 4],
     [8, 9, 3, 1, 0, 6, 4, 2, 5, 7],
     [1, 7, 6, 5, 4, 8, 0, 3, 2, 9],
     [2, 1, 9, 0, 6, 7, 5, 8, 4, 3],
     [5, 2, 7, 4, 9, 3, 1, 0, 8, 6],
     [4, 3, 0, 6, 1, 5, 2, 9, 7, 8],
     [9, 8, 5, 7, 2, 0, 3, 4, 6, 1],
     [7, 0, 1, 9, 5, 4, 8, 6, 3, 2],
     [6, 5, 4, 3, 8, 2, 7, 1, 9, 0]]


The method ``is_division_algebra`` tests for the cancellation property,
as shown below.

.. code:: ipython3

    >>> latin_sqr.has_cancellation()




.. parsed-literal::

    True



.. code:: ipython3

    >>> latin_sqr.is_associative()




.. parsed-literal::

    False



Loops
~~~~~

A **Loop** is a Quasigroup with an identity element.

The following example is from *“Counting loops with the inverse
property”* by Asif Ali & John Slaney
(http://www.quasigroups.eu/contents/download/2008/16_2.pdf)

.. code:: ipython3

    >>> ip_loop = fa.make_finite_algebra('IP_Loop',
                                         'IP loop of small order',
                                         ('0', '1', '2', '3', '4', '5', '6'),
                                         [[0, 1, 2, 3, 4, 5, 6], 
                                          [1, 2, 0, 5, 6, 4, 3], 
                                          [2, 0, 1, 6, 5, 3, 4], 
                                          [3, 6, 5, 4, 0, 1, 2], 
                                          [4, 5, 6, 0, 3, 2, 1], 
                                          [5, 3, 4, 2, 1, 6, 0], 
                                          [6, 4, 3, 1, 2, 0, 5]])
    >>> ip_loop.about()


.. parsed-literal::

    
    ** Loop **
    Name: IP_Loop
    Instance ID: 5215432976
    Description: IP loop of small order
    Order: 7
    Identity: 0
    Associative? No
    Commutative? No
    Cyclic?: No
    Elements: ('0', '1', '2', '3', '4', '5', '6')
    Has Cancellation? Yes
    Has Inverses? Yes
    Cayley Table (showing indices):
    [[0, 1, 2, 3, 4, 5, 6],
     [1, 2, 0, 5, 6, 4, 3],
     [2, 0, 1, 6, 5, 3, 4],
     [3, 6, 5, 4, 0, 1, 2],
     [4, 5, 6, 0, 3, 2, 1],
     [5, 3, 4, 2, 1, 6, 0],
     [6, 4, 3, 1, 2, 0, 5]]


.. code:: ipython3

    >>> ip_loop.has_cancellation()




.. parsed-literal::

    True



.. code:: ipython3

    >>> ip_loop.identity




.. parsed-literal::

    '0'



.. code:: ipython3

    >>> ip_loop.is_associative()




.. parsed-literal::

    False



Semigroup
~~~~~~~~~

**Semigroup** – an associative Magma: for any
:math:`a,b,c \in S \Rightarrow a \circ (b \circ c) = (a \circ b) \circ c`

**Example:**

Reference: `Groupoids and Smarandache
Groupoids <https://arxiv.org/ftp/math/papers/0304/0304490.pdf>`__ by W.
B. Vasantha Kandasamy

.. code:: ipython3

    >>> sg = fa.make_finite_algebra(
        'Example 1.4.1',
        'See: Groupoids and Smarandache Groupoids by W. B. Vasantha Kandasamy',
        ['a', 'b', 'c', 'd', 'e', 'f'],
        [[0, 3, 0, 3, 0, 3],
         [1, 4, 1, 4, 1, 4],
         [2, 5, 2, 5, 2, 5],
         [3, 0, 3, 0, 3, 0],
         [4, 1, 4, 1, 4, 1],
         [5, 2, 5, 2, 5, 2]]
    )
    
    >>> sg.about()


.. parsed-literal::

    
    ** Semigroup **
    Name: Example 1.4.1
    Instance ID: 5215434896
    Description: See: Groupoids and Smarandache Groupoids by W. B. Vasantha Kandasamy
    Order: 6
    Identity: None
    Associative? Yes
    Commutative? No
    Cyclic?: No
    Elements: ('a', 'b', 'c', 'd', 'e', 'f')
    Has Cancellation? No
    Has Inverses? No
    Cayley Table (showing indices):
    [[0, 3, 0, 3, 0, 3],
     [1, 4, 1, 4, 1, 4],
     [2, 5, 2, 5, 2, 5],
     [3, 0, 3, 0, 3, 0],
     [4, 1, 4, 1, 4, 1],
     [5, 2, 5, 2, 5, 2]]


Since the element in the 0,1 position of the table is 3, it follows
that, :math:`a \circ b = d`:

.. code:: ipython3

    >>> sg.op('a', 'b')




.. parsed-literal::

    'd'



This Semigroup is regular and every element has three weak inverses.

.. code:: ipython3

    >>> sg.is_regular()




.. parsed-literal::

    True



.. code:: ipython3

    >>> sg.weak_inverses()




.. parsed-literal::

    {'a': ['a', 'c', 'e'],
     'b': ['b', 'd', 'f'],
     'c': ['a', 'c', 'e'],
     'd': ['b', 'd', 'f'],
     'e': ['a', 'c', 'e'],
     'f': ['b', 'd', 'f']}



By the way, the output of the ``about`` method, above, indicated that
there are 5 more **generators** in addition to the 2 that were printed
out. To obtain the full list of generators simply apply the
``generators`` method to the algebra, as shown below. In this case, the
entire algebra can be obtained by computing the closure of any of one of
the collections of 3 generators.

.. code:: ipython3

    >>> help(sg.generators)


.. parsed-literal::

    Help on method generators in module finalg.magma:
    
    generators(start_of_range=1) method of finalg.semigroup.Semigroup instance
        If the algebra is cyclic, then a list of individual elements that each
        generate the algebra is returned; otherwise, a list of lists of elements,
        is returned, where each sublist generates the algebra. This method looks
        for the smallest sets of elements that can generate the group. It stops
        looking once it finds all small sets of elements of a given size.
    


.. code:: ipython3

    >>> gens = sg.generators()
    >>> print("\nGenerators:")
    >>> gens


.. parsed-literal::

    
    Generators:




.. parsed-literal::

    [('a', 'b', 'c'),
     ('a', 'b', 'f'),
     ('a', 'e', 'f'),
     ('b', 'c', 'd'),
     ('b', 'd', 'f'),
     ('c', 'd', 'e'),
     ('d', 'e', 'f')]



.. code:: ipython3

    >>> for gen in gens:
    >>>     print(f"sg.closure{gen} = {sg.closure(gen, include_inverses=False)}")


.. parsed-literal::

    sg.closure('a', 'b', 'c') = ['e', 'a', 'b', 'f', 'd', 'c']
    sg.closure('a', 'b', 'f') = ['e', 'a', 'b', 'f', 'd', 'c']
    sg.closure('a', 'e', 'f') = ['e', 'a', 'b', 'f', 'd', 'c']
    sg.closure('b', 'c', 'd') = ['e', 'a', 'b', 'f', 'd', 'c']
    sg.closure('b', 'd', 'f') = ['e', 'a', 'b', 'f', 'd', 'c']
    sg.closure('c', 'd', 'e') = ['e', 'a', 'b', 'f', 'd', 'c']
    sg.closure('d', 'e', 'f') = ['e', 'a', 'b', 'f', 'd', 'c']


Monoid
~~~~~~

**Monoid** – a Semigroup with identity element: :math:`\exists e \in S`,
such that, for all :math:`a \in S, a \circ e = e \circ a = a`

.. code:: ipython3

    >>> m4 = fa.make_finite_algebra(
        'M4',
        'Example of a commutative monoid',
        ['a', 'b', 'c', 'd'],
        [[0, 0, 0, 0],
         [0, 1, 2, 3],
         [0, 2, 0, 2],
         [0, 3, 2, 1]])
    
    >>> m4.about(use_table_names=True)


.. parsed-literal::

    
    ** Monoid **
    Name: M4
    Instance ID: 5215432656
    Description: Example of a commutative monoid
    Order: 4
    Identity: b
    Associative? Yes
    Commutative? Yes
    Cyclic?: No
    Elements: ('a', 'b', 'c', 'd')
    Has Cancellation? No
    Has Inverses? No
    Cayley Table (showing names):
    [['a', 'a', 'a', 'a'],
     ['a', 'b', 'c', 'd'],
     ['a', 'c', 'a', 'c'],
     ['a', 'd', 'c', 'b']]


By the way, the Monoid, above, and others like it of different orders,
can be automatically generated using the function,
``generate_commutative_monoid``. It is based on integer multiplication
modulo the desired order.

.. code:: ipython3

    >>> m4.identity  # Returns the identity element




.. parsed-literal::

    'b'



.. code:: ipython3

    >>> m4.op('c', 'b')  # since 'b' is the identity element




.. parsed-literal::

    'c'



Serialization
-------------

Algebras can be converted to and from JSON strings/files and Python
dictionaries.

Instantiate Algebra from JSON File
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

First setup some path variables:

-  one that points to the abstract_algebra directory
-  and the other points to a subdirectory containing algebra definitions
   in JSON format

Also, the code here assumes that there is an environment variable,
``PYPROJ``, that points to the parent directory of the abstract_algebra
directory.

.. code:: ipython3

    >>> import os
    >>> aa_path = os.path.join(os.getenv("PYPROJ"), "abstract_algebra")
    >>> alg_dir = os.path.join(aa_path, "Algebras")

Here’s the **JSON file**:

The path to the JSON file is constructed by using Python’s
*os.path.join* to join strings together. So, the quantity, ``v4_json``,
below, is a string. And then we’ve used Jupyter Notebook’s ability to
“reach out” into the external environment (via “!”) and printout the
file using the UNIX command, ``cat``.

.. code:: ipython3

    >>> v4_json = os.path.join(alg_dir, "v4_klein_4_group.json")
    
    >>> !cat {v4_json}


.. parsed-literal::

    {"name": "V4",
     "description": "Klein-4 group",
     "elements": ["e", "h", "v", "r"],
     "table": [[0, 1, 2, 3],
               [1, 0, 3, 2],
               [2, 3, 0, 1],
               [3, 2, 1, 0]]
    }


And, here’s the **algebra** that is loaded from the JSON file:

.. code:: ipython3

    >>> v4 = fa.make_finite_algebra(v4_json)
    
    >>> v4




.. parsed-literal::

    Group(
    'V4',
    'Klein-4 group',
    ('e', 'h', 'v', 'r'),
    [[0, 1, 2, 3], [1, 0, 3, 2], [2, 3, 0, 1], [3, 2, 1, 0]]
    )



Convert Algebra to Python Dictionary
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The examples, below, show a Magma and a Group being converted into
dictionaries.

.. code:: ipython3

    >>> rps.to_dict()




.. parsed-literal::

    {'name': 'RPS',
     'description': 'Rock, Paper, Scissors Magma',
     'elements': ('r', 'p', 's'),
     'table': [[0, 1, 0], [1, 1, 2], [0, 2, 2]]}



Note that the only difference between the JSON and Python dictionary
representations of an algebra is the type of quotes used aroung strings.
JSON requires that double quotes be used, while Python uses single
quotes by default.

In the second example, below, the **type** of algebra (e.g., Magma) can
be included in the dictionary for readability, however, the *type* field
is ignored when ``make_finite_algebra`` reads a dictionary or JSON file.
Again, the type of algebra and its properties are always derived from
its Cayley table.

.. code:: ipython3

    >>> rps_dict = rps.to_dict(include_classname=True)
    
    >>> rps_dict




.. parsed-literal::

    {'name': 'RPS',
     'description': 'Rock, Paper, Scissors Magma',
     'elements': ('r', 'p', 's'),
     'table': [[0, 1, 0], [1, 1, 2], [0, 2, 2]],
     'type': 'Magma'}



.. code:: ipython3

    >>> v4_dict = v4.to_dict()
    
    >>> v4_dict




.. parsed-literal::

    {'name': 'V4',
     'description': 'Klein-4 group',
     'elements': ('e', 'h', 'v', 'r'),
     'table': [[0, 1, 2, 3], [1, 0, 3, 2], [2, 3, 0, 1], [3, 2, 1, 0]]}



Instantiate Algebra from Python Dictionary
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

For these examples, we’ll use the Python dictionaries, created above, as
inputs to ``make_finite_algebra``.

.. code:: ipython3

    >>> rps_from_dict = fa.make_finite_algebra(rps_dict)
    
    >>> rps_from_dict




.. parsed-literal::

    Magma(
    'RPS',
    'Rock, Paper, Scissors Magma',
    ('r', 'p', 's'),
    [[0, 1, 0], [1, 1, 2], [0, 2, 2]]
    )



.. code:: ipython3

    >>> v4_from_dict = fa.make_finite_algebra(v4_dict)
    
    >>> v4_from_dict




.. parsed-literal::

    Group(
    'V4',
    'Klein-4 group',
    ('e', 'h', 'v', 'r'),
    [[0, 1, 2, 3], [1, 0, 3, 2], [2, 3, 0, 1], [3, 2, 1, 0]]
    )



Convert Algebra to JSON String
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Note that the conversion example here outputs a single Python string
(i.e., enclosed by single quotes), within which, all the strings are
enclosed by double quotes, as required by JSON.

.. code:: ipython3

    >>> v4_json_string = v4.dumps()
    
    >>> v4_json_string




.. parsed-literal::

    '{"name": "V4", "description": "Klein-4 group", "elements": ["e", "h", "v", "r"], "table": [[0, 1, 2, 3], [1, 0, 3, 2], [2, 3, 0, 1], [3, 2, 1, 0]]}'



**WARNING**: Although an algebra can be constructed by loading its
definition from a JSON file, it cannot be constructed directly from a
JSON string, because ``make_finite_algebra`` interprets a single string
input as a JSON file name. To load an algebra from a JSON string, first
convert the string to a Python dictionary, then input that to
``make_finite_algebra``, as shown below, using the JSON string
constructed above:

.. code:: ipython3

    >>> import json
    
    >>> fa.make_finite_algebra(json.loads(v4_json_string))




.. parsed-literal::

    Group(
    'V4',
    'Klein-4 group',
    ('e', 'h', 'v', 'r'),
    [[0, 1, 2, 3], [1, 0, 3, 2], [2, 3, 0, 1], [3, 2, 1, 0]]
    )



Autogeneration of Finite Algebras
---------------------------------

There are several functions for autogenerating finite algebras of any
desired size:

**Groups**

-  ``generate_cyclic_group(n)``: :math:`Z_n`, where
   :math:`a \circ b \equiv a+b` mod :math:`n`, where
   :math:`a,b \in \{0,1,...,n-1\}`; order is :math:`n`
-  ``generate_symmetric_group(n)``: :math:`S_n`, where :math:`\circ` is
   composition of permutations of :math:`(0, 1, ..., n-1)`; order is
   :math:`n!`
-  ``generate_powerset_group(n)``:
   :math:`A \circ B \equiv A \bigtriangleup B`, where
   :math:`A,B \in P(\{0, 1, ..., n-1\})`; order is :math:`2^n`

**Monoid**

-  ``generate_commutative_monoid(n)``: :math:`a \circ b \equiv ab` mod
   :math:`n`, where :math:`a,b \in \{0,1,...,n-1\}`; order is :math:`n`

Autogenerated Cyclic Group
~~~~~~~~~~~~~~~~~~~~~~~~~~

A cyclic group of any desired order can be generated. A very small one
is created, below, because it will be used later to demonstrate Direct
Products and Isomorphisms.

.. code:: ipython3

    >>> z2 = fa.generate_cyclic_group(2)
    >>> z2.about()


.. parsed-literal::

    
    ** Group **
    Name: Z2
    Instance ID: 5216539888
    Description: Autogenerated cyclic Group of order 2
    Order: 2
    Identity: '0'
    Commutative? Yes
    Cyclic?: Yes
    Generators: ['1']
    Elements:
       Index   Name   Inverse  Order
          0     '0'     '0'       1
          1     '1'     '1'       2
    Cayley Table (showing indices):
    [[0, 1], [1, 0]]




.. parsed-literal::

    '<Group:Z2, ID:5216539888>'



Autogenerated Symmetric Group
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The symmetric group, based on the permutations of n elements, (1, 2, 3,
…, n), can be generated as follows:

**WARNING**: Since the order of an autogenerated symmetric group is
**n!**, even a small value of **n** can result in a very large group.

.. code:: ipython3

    >>> s3 = fa.generate_symmetric_group(3)
    >>> s3.about()


.. parsed-literal::

    
    ** Group **
    Name: S3
    Instance ID: 5216451440
    Description: Autogenerated symmetric Group on 3 elements
    Order: 6
    Identity: '(0, 1, 2)'
    Commutative? No
    Cyclic?: No
    Elements:
       Index   Name   Inverse  Order
          0 '(0, 1, 2)' '(0, 1, 2)'       1
          1 '(0, 2, 1)' '(0, 2, 1)'       2
          2 '(1, 0, 2)' '(1, 0, 2)'       2
          3 '(1, 2, 0)' '(2, 0, 1)'       3
          4 '(2, 0, 1)' '(1, 2, 0)'       3
          5 '(2, 1, 0)' '(2, 1, 0)'       2
    Cayley Table (showing indices):
    [[0, 1, 2, 3, 4, 5],
     [1, 0, 4, 5, 2, 3],
     [2, 3, 0, 1, 5, 4],
     [3, 2, 5, 4, 0, 1],
     [4, 5, 1, 0, 3, 2],
     [5, 4, 3, 2, 1, 0]]




.. parsed-literal::

    '<Group:S3, ID:5216451440>'



Autogenerated Powerset Group
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The function, ``generate_powerset_group``, will generate a group on the
powerset of {0, 1, 2, …, n-1} with **symmetric difference** as the
group’s binary operation. This group is useful because it can be used to
form a ring with set intersection as the second operator.

This means that the order of the autogenerated powerset group will be
:math:`2^n`, so the same WARNING as above applies with regard to large
values of n.

.. code:: ipython3

    >>> ps3 = fa.generate_powerset_group(3)
    >>> ps3.about()


.. parsed-literal::

    
    ** Group **
    Name: PS3
    Instance ID: 5216356784
    Description: Autogenerated Group on the powerset of 3 elements, with symmetric difference operator
    Order: 8
    Identity: '{}'
    Commutative? Yes
    Cyclic?: No
    Elements:
       Index   Name   Inverse  Order
          0    '{}'    '{}'       1
          1   '{0}'   '{0}'       2
          2   '{1}'   '{1}'       2
          3   '{2}'   '{2}'       2
          4 '{0, 1}' '{0, 1}'       2
          5 '{0, 2}' '{0, 2}'       2
          6 '{1, 2}' '{1, 2}'       2
          7 '{0, 1, 2}' '{0, 1, 2}'       2
    Cayley Table (showing indices):
    [[0, 1, 2, 3, 4, 5, 6, 7],
     [1, 0, 4, 5, 2, 3, 7, 6],
     [2, 4, 0, 6, 1, 7, 3, 5],
     [3, 5, 6, 0, 7, 1, 2, 4],
     [4, 2, 1, 7, 0, 6, 5, 3],
     [5, 3, 7, 1, 6, 0, 4, 2],
     [6, 7, 3, 2, 5, 4, 0, 1],
     [7, 6, 5, 4, 3, 2, 1, 0]]




.. parsed-literal::

    '<Group:PS3, ID:5216356784>'



Autogenerated Monoid
~~~~~~~~~~~~~~~~~~~~

The function, ``generate_commutative_monoid``, is based on integer
multiplication modulo the desired order.

.. code:: ipython3

    >>> m7 = fa.generate_commutative_monoid(7)
    >>> m7.about()


.. parsed-literal::

    
    ** Monoid **
    Name: M7
    Instance ID: 5215724528
    Description: Autogenerated commutative Monoid of order 7
    Order: 7
    Identity: a1
    Associative? Yes
    Commutative? Yes
    Cyclic?: No
    Elements: ('a0', 'a1', 'a2', 'a3', 'a4', 'a5', 'a6')
    Has Cancellation? No
    Has Inverses? No
    Cayley Table (showing indices):
    [[0, 0, 0, 0, 0, 0, 0],
     [0, 1, 2, 3, 4, 5, 6],
     [0, 2, 4, 6, 1, 3, 5],
     [0, 3, 6, 2, 5, 1, 4],
     [0, 4, 1, 5, 2, 6, 3],
     [0, 5, 3, 1, 6, 4, 2],
     [0, 6, 5, 4, 3, 2, 1]]


Direct Products
---------------

The **direct product** of two or more algebras can be generated using
Python’s multiplication operator, ``*``:

Direct Product of Multiple Groups
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code:: ipython3

    >>> z2_cubed = z2 * z2 * z2
    
    >>> z2_cubed.about()


.. parsed-literal::

    
    ** Group **
    Name: Z2_x_Z2_x_Z2
    Instance ID: 5182937424
    Description: Direct product of Z2_x_Z2 & Z2
    Order: 8
    Identity: '0:0:0'
    Commutative? Yes
    Cyclic?: No
    Elements:
       Index   Name   Inverse  Order
          0 '0:0:0' '0:0:0'       1
          1 '0:0:1' '0:0:1'       2
          2 '0:1:0' '0:1:0'       2
          3 '0:1:1' '0:1:1'       2
          4 '1:0:0' '1:0:0'       2
          5 '1:0:1' '1:0:1'       2
          6 '1:1:0' '1:1:0'       2
          7 '1:1:1' '1:1:1'       2
    Cayley Table (showing indices):
    [[0, 1, 2, 3, 4, 5, 6, 7],
     [1, 0, 3, 2, 5, 4, 7, 6],
     [2, 3, 0, 1, 6, 7, 4, 5],
     [3, 2, 1, 0, 7, 6, 5, 4],
     [4, 5, 6, 7, 0, 1, 2, 3],
     [5, 4, 7, 6, 1, 0, 3, 2],
     [6, 7, 4, 5, 2, 3, 0, 1],
     [7, 6, 5, 4, 3, 2, 1, 0]]




.. parsed-literal::

    '<Group:Z2_x_Z2_x_Z2, ID:5182937424>'



Direct Product of Monoids
~~~~~~~~~~~~~~~~~~~~~~~~~

.. code:: ipython3

    >>> mon3 = fa.generate_commutative_monoid(3)
    >>> mon3




.. parsed-literal::

    Monoid(
    'M3',
    'Autogenerated commutative Monoid of order 3',
    ('a0', 'a1', 'a2'),
    [[0, 0, 0], [0, 1, 2], [0, 2, 1]]
    )



.. code:: ipython3

    >>> m3_sqr = mon3 * mon3
    >>> m3_sqr.about()


.. parsed-literal::

    
    ** Monoid **
    Name: M3_x_M3
    Instance ID: 5215417936
    Description: Direct product of M3 & M3
    Order: 9
    Identity: a1:a1
    Associative? Yes
    Commutative? Yes
    Cyclic?: No
    Elements: ('a0:a0', 'a0:a1', 'a0:a2', 'a1:a0', 'a1:a1', 'a1:a2', 'a2:a0', 'a2:a1', 'a2:a2')
    Has Cancellation? No
    Has Inverses? No
    Cayley Table (showing indices):
    [[0, 0, 0, 0, 0, 0, 0, 0, 0],
     [0, 1, 2, 0, 1, 2, 0, 1, 2],
     [0, 2, 1, 0, 2, 1, 0, 2, 1],
     [0, 0, 0, 3, 3, 3, 6, 6, 6],
     [0, 1, 2, 3, 4, 5, 6, 7, 8],
     [0, 2, 1, 3, 5, 4, 6, 8, 7],
     [0, 0, 0, 6, 6, 6, 3, 3, 3],
     [0, 1, 2, 6, 7, 8, 3, 4, 5],
     [0, 2, 1, 6, 8, 7, 3, 5, 4]]


Isomorphisms
------------

If two algebras are isomorphic, then the mapping between their elements
is returned as a Python dictionary.

Here is a well-known example, using two small groups created above, v4
and the direct product of z2 with itself, z2 \* z2:

Group Isomorphism
~~~~~~~~~~~~~~~~~

.. code:: ipython3

    >>> z2_sqr = z2 * z2
    
    >>> v4.isomorphic(z2_sqr)




.. parsed-literal::

    {'e': '0:0', 'h': '0:1', 'v': '1:0', 'r': '1:1'}



If two algebras are not isomorphic, then ``False`` is returned.

.. code:: ipython3

    >>> z4 = fa.generate_cyclic_group(4)
    
    >>> z4.isomorphic(z2_sqr)




.. parsed-literal::

    False



Magma Isomorphism
~~~~~~~~~~~~~~~~~

In this example, we’ll use a made-up Magma, similar to Rock, Paper,
Scissors.

**Water, Fire, Stick:**

-  Water quenches Fire
-  Fire burns Stick
-  Stick floats on Water

.. code:: ipython3

    >>> wfs = fa.make_finite_algebra(
        'WFS',
        'Water, Fire, Stick Magma',
        ['water', 'fire', 'stick'],
        [[0, 0, 2],
         [0, 1, 1],
         [2, 1, 2]])
    
    >>> wfs




.. parsed-literal::

    Magma(
    'WFS',
    'Water, Fire, Stick Magma',
    ('water', 'fire', 'stick'),
    [[0, 0, 2], [0, 1, 1], [2, 1, 2]]
    )



Here’s the isomorphism between rps and wfs:

.. code:: ipython3

    >>> rps.isomorphic(wfs)




.. parsed-literal::

    {'r': 'water', 'p': 'stick', 's': 'fire'}



Subalgebras (Subgroups)
-----------------------

An algebra can contain subalgebras (e.g., a Group can have subgroups).
In fact, sometimes the subalgebra may not even be of the same type as
the parent algebra. For, example, we’ll see below that a Semigroup can
contain a Group as a subalgebra.

The method, ``proper_subalgebras``, extracts all possible proper
subalgebras that exist within an algebra, regardless of whether they are
isomorphic to each other or not, or even of the same algebraic class as
the parent algebra.

A subalgebra is constructed from a subset of elements that are closed
under the algebra’s binary operation.

Example: Proper Subgroups
~~~~~~~~~~~~~~~~~~~~~~~~~

.. code:: ipython3

    >>> z8 = fa.generate_cyclic_group(8)
    >>> z8.about()


.. parsed-literal::

    
    ** Group **
    Name: Z8
    Instance ID: 5216939024
    Description: Autogenerated cyclic Group of order 8
    Order: 8
    Identity: '0'
    Commutative? Yes
    Cyclic?: Yes
    Generators: ['7', '1', '3', '5']
    Elements:
       Index   Name   Inverse  Order
          0     '0'     '0'       1
          1     '1'     '7'       8
          2     '2'     '6'       4
          3     '3'     '5'       8
          4     '4'     '4'       2
          5     '5'     '3'       8
          6     '6'     '2'       4
          7     '7'     '1'       8
    Cayley Table (showing indices):
    [[0, 1, 2, 3, 4, 5, 6, 7],
     [1, 2, 3, 4, 5, 6, 7, 0],
     [2, 3, 4, 5, 6, 7, 0, 1],
     [3, 4, 5, 6, 7, 0, 1, 2],
     [4, 5, 6, 7, 0, 1, 2, 3],
     [5, 6, 7, 0, 1, 2, 3, 4],
     [6, 7, 0, 1, 2, 3, 4, 5],
     [7, 0, 1, 2, 3, 4, 5, 6]]




.. parsed-literal::

    '<Group:Z8, ID:5216939024>'



.. code:: ipython3

    >>> z8_proper_subs = z8.proper_subalgebras()
    
    >>> for sub in z8_proper_subs:
    >>>     sub.about()


.. parsed-literal::

    
    ** Group **
    Name: Z8_subalgebra_0
    Instance ID: 5216863568
    Description: Subalgebra of: Autogenerated cyclic Group of order 8
    Order: 2
    Identity: '0'
    Commutative? Yes
    Cyclic?: Yes
    Generators: ['4']
    Elements:
       Index   Name   Inverse  Order
          0     '0'     '0'       1
          1     '4'     '4'       2
    Cayley Table (showing indices):
    [[0, 1], [1, 0]]
    
    ** Group **
    Name: Z8_subalgebra_1
    Instance ID: 5216863248
    Description: Subalgebra of: Autogenerated cyclic Group of order 8
    Order: 4
    Identity: '0'
    Commutative? Yes
    Cyclic?: Yes
    Generators: ['6', '2']
    Elements:
       Index   Name   Inverse  Order
          0     '0'     '0'       1
          1     '2'     '6'       4
          2     '4'     '4'       2
          3     '6'     '2'       4
    Cayley Table (showing indices):
    [[0, 1, 2, 3], [1, 2, 3, 0], [2, 3, 0, 1], [3, 0, 1, 2]]


Normal Subgroups
~~~~~~~~~~~~~~~~

Both of the subgroups of Z8, derived above, are **normal**:

.. code:: ipython3

    >>> [z8.is_normal(g) for g in z8_proper_subs]




.. parsed-literal::

    [True, True]



Proper Subalgebras up to Isomorphism
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The function, ``about_subalgebras``, finds all the subalgebras (e.g.,
subgroups) of an algebra, and then partitions that list into sublists
(referred to as partitions here) where the algebras in each partition
are isomporphic to each other. It returns the list of partitions and
prints a summary of each partition.

The example, below, demonstrats this using the autogenerated powerset
group, **ps3**, that was created earlier.

.. code:: ipython3

    >>> iso_parts = fa.about_subalgebras(ps3)


.. parsed-literal::

    
    Subalgebras of PS3 : Autogenerated Group on the powerset of 3 elements, with symmetric difference operator
    
      There are 2 unique proper subalgebras, up to isomorphism, out of 14 total subalgebras.
      as shown by the partitions below:
    
    7 Isomorphic Commutative Normal Groups of order 4 with identity '{}':
          Group: PS3_subalgebra_0: ('{}', '{0, 1}', '{0, 2}', '{1, 2}')
          Group: PS3_subalgebra_3: ('{}', '{2}', '{0, 1}', '{0, 1, 2}')
          Group: PS3_subalgebra_5: ('{}', '{0}', '{1}', '{0, 1}')
          Group: PS3_subalgebra_6: ('{}', '{0}', '{2}', '{0, 2}')
          Group: PS3_subalgebra_8: ('{}', '{1}', '{0, 2}', '{0, 1, 2}')
          Group: PS3_subalgebra_9: ('{}', '{0}', '{1, 2}', '{0, 1, 2}')
          Group: PS3_subalgebra_11: ('{}', '{1}', '{2}', '{1, 2}')
    
    7 Isomorphic Commutative Normal Groups of order 2 with identity '{}':
          Group: PS3_subalgebra_1: ('{}', '{0}')
          Group: PS3_subalgebra_2: ('{}', '{0, 1}')
          Group: PS3_subalgebra_4: ('{}', '{1}')
          Group: PS3_subalgebra_7: ('{}', '{1, 2}')
          Group: PS3_subalgebra_10: ('{}', '{0, 1, 2}')
          Group: PS3_subalgebra_12: ('{}', '{0, 2}')
          Group: PS3_subalgebra_13: ('{}', '{2}')
    


And here, for example, is the first subalgebra found in the first
partition:

.. code:: ipython3

    >>> iso_parts[0][0]




.. parsed-literal::

    Group(
    'PS3_subalgebra_0',
    'Subalgebra of: Autogenerated Group on the powerset of 3 elements, with symmetric difference operator',
    ('{}', '{0, 1}', '{0, 2}', '{1, 2}'),
    [[0, 1, 2, 3], [1, 0, 3, 2], [2, 3, 0, 1], [3, 2, 1, 0]]
    )



Subalgebras of Semigroups, Etc.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Recall the Semigroup example from above:

.. code:: ipython3

    >>> sg.about()


.. parsed-literal::

    
    ** Semigroup **
    Name: Example 1.4.1
    Instance ID: 5215434896
    Description: See: Groupoids and Smarandache Groupoids by W. B. Vasantha Kandasamy
    Order: 6
    Identity: None
    Associative? Yes
    Commutative? No
    Cyclic?: No
    Elements: ('a', 'b', 'c', 'd', 'e', 'f')
    Has Cancellation? No
    Has Inverses? No
    Cayley Table (showing indices):
    [[0, 3, 0, 3, 0, 3],
     [1, 4, 1, 4, 1, 4],
     [2, 5, 2, 5, 2, 5],
     [3, 0, 3, 0, 3, 0],
     [4, 1, 4, 1, 4, 1],
     [5, 2, 5, 2, 5, 2]]


As we will see, below, the Semigroup, sg, contains 4 unique subalgebras,
up to isomorphism:

-  3 Semigroups and
-  1 Group

.. code:: ipython3

    >>> iso_parts = fa.about_subalgebras(sg)


.. parsed-literal::

    
    Subalgebras of Example 1.4.1 : See: Groupoids and Smarandache Groupoids by W. B. Vasantha Kandasamy
    
      There are 4 unique proper subalgebras, up to isomorphism, out of 10 total subalgebras.
      as shown by the partitions below:
    
    3 Isomorphic Semigroups of order 4:
          Semigroup: Example 1.4.1_subalgebra_0: ('a', 'b', 'd', 'e')
          Semigroup: Example 1.4.1_subalgebra_3: ('a', 'c', 'd', 'f')
          Semigroup: Example 1.4.1_subalgebra_6: ('b', 'c', 'e', 'f')
    
    3 Isomorphic Commutative Groups of order 2:
          Group: Example 1.4.1_subalgebra_1: ('a', 'd') with identity 'a'
          Group: Example 1.4.1_subalgebra_2: ('c', 'f') with identity 'c'
          Group: Example 1.4.1_subalgebra_8: ('b', 'e') with identity 'e'
    
    3 Isomorphic Semigroups of order 2:
          Semigroup: Example 1.4.1_subalgebra_4: ('c', 'e')
          Semigroup: Example 1.4.1_subalgebra_7: ('a', 'e')
          Semigroup: Example 1.4.1_subalgebra_9: ('a', 'c')
    
    1 Semigroup of order 3:
          Semigroup: Example 1.4.1_subalgebra_5: ('a', 'c', 'e')
    

