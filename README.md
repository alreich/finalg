# Finite Abstract Algebras (finalg)

Finite abstract algebras in Python - groups, rings, fields, vector spaces, etc.

*(This is a continuation of the now archived repo, https://github.com/alreich/abstract_algebra)*

**finalg** is a Python package for constructing and exploring **finite abstract algebras** —
Magmas, Semigroups, Monoids, Groups, Quasigroups, Loops, Rings, Fields, Modules, and Vector
Spaces — using explicit **Cayley (multiplication) tables**. It is aimed at hands-on exploration
of abstract algebra: building small algebras by hand or by generator function, inspecting their
properties (associativity, commutativity, identities, inverses, generators, subalgebras, ...),
combining them (direct products, quotient groups, Cayley–Dickson construction), and representing
them in other forms (regular representations, abstract matrices, Cayley diagrams).

- Source: <https://github.com/alreich/finalg>
- Documentation: <https://finalg.readthedocs.io/en/latest/>

## Table of Contents

- [Core Ideas](#core-ideas)
- [The Algebra Class Hierarchy](#the-algebra-class-hierarchy)
- [Elements, Cayley Tables, and Operators](#elements-cayley-tables-and-operators)
- [Creating Algebras](#creating-algebras)
- [Serialization](#serialization)
- [Arithmetic on Elements](#arithmetic-on-elements)
- [Algebra-Level Operators: Direct Products and Quotients](#algebra-level-operators-direct-products-and-quotients)
- [What Each Level of the Hierarchy Adds](#what-each-level-of-the-hierarchy-adds)
- [Rings and Fields](#rings-and-fields)
- [The Cayley–Dickson Construction](#the-cayley-dickson-construction)
- [Modules and Vector Spaces](#modules-and-vector-spaces)
- [Automatically Generating Algebras](#automatically-generating-algebras)
- [Supporting Modules](#supporting-modules)
- [Repository Layout](#repository-layout)
- [References](#references)

## Core Ideas

- Every algebra's elements are **strings**. This keeps algebras easy to read, print, and
  serialize, and lets element names carry meaning (e.g., permutations in cycle notation,
  cosets prefixed with `~`, direct-product elements joined with a delimiter like `e:a`).
- Every algebra with a single binary operation is backed by a **Cayley table**: an NxN table
  (`N` = number of elements) where entry `(row, col)` gives the result of applying the operator
  to the row element and the column element, in that order (`row * col`). Rings and Fields carry
  a second Cayley table for their second operator.
- Properties like associativity, commutativity, identity, and inverses are **derived entirely
  from the table(s)**, not asserted by the class. Creating an algebra with `make_finite_algebra`
  inspects the table and returns the most specific class the table supports.
- All of an algebra's structural and single-element-set classes ultimately share one abstract
  base, `FiniteAlgebra`. Multi-element-set algebras (Modules, Vector Spaces) share a separate
  base, `FiniteCompositeAlgebra`.

## The Algebra Class Hierarchy

```
FiniteAlgebra (ABC)
└── Magma
    ├── Semigroup
    │   └── Monoid
    │       └── Group
    │           └── Ring
    │               └── Field
    └── Quasigroup
        └── Loop

FiniteCompositeAlgebra (ABC)
└── Module
    └── VectorSpace

NDimensionalModule(Module)         # constructed directly from a Ring
NDimensionalVectorSpace(VectorSpace)  # constructed directly from a Field
```

`FiniteAlgebra` and `FiniteCompositeAlgebra` are not intended to be instantiated directly; the
highest classes in each hierarchy meant for direct use are `Magma` and `Module`/`VectorSpace`,
respectively.

The general design strategy is to place each method **as high in the hierarchy as its
requirements allow**. `Magma` only requires a binary operator, so most generic structural
methods (direct products, closures, subalgebra discovery, isomorphism checking, Cayley
diagrams, ...) live there. `Semigroup` adds associativity, so associativity-dependent methods
live no higher than `Semigroup`. `Monoid` adds an identity element, so identity-dependent
methods (e.g., element order) live no higher than `Monoid`. `Group` adds inverses, so
inverse-dependent methods (conjugation, commutators, cosets, quotient groups, solvability) live
no higher than `Group`. `Ring` adds a second, distributive operator, and `Field` adds that
every non-zero element has a multiplicative inverse.

`Quasigroup` and `Loop` branch off of `Magma` independently of the `Semigroup`/`Monoid`/`Group`
line: a `Quasigroup` supports unique left/right division (equivalently, cancellation) without
necessarily being associative, and a `Loop` is a `Quasigroup` that also has an identity element.

Constructors for the individual classes exist but are **not meant to be called directly**.
Instead, use the single factory function `make_finite_algebra` (see
[Creating Algebras](#creating-algebras)), which inspects the input table(s) and returns an
instance of the most specific applicable class:

| Table 1 property | Table 2? | Result |
|---|---|---|
| Not associative, no cancellation | — | `Magma` |
| Not associative, cancellation, no identity | — | `Quasigroup` |
| Not associative, cancellation, identity | — | `Loop` |
| Associative, no identity | — | `Semigroup` |
| Associative, identity, no inverses | — | `Monoid` |
| Associative, identity, inverses | No | `Group` |
| Associative, identity, inverses | Yes, associative & distributive | `Ring` (or `Field`, if elements minus the additive identity form a commutative multiplicative group) |

## Elements, Cayley Tables, and Operators

- **`CayleyTable`** (`cayley_table.py`) wraps a square table of element indices (a NumPy array
  under the hood) and provides the property checks described above (`is_associative`,
  `is_commutative`, `has_cancellation`, `has_inverses`, `identity`, `distributes_over`, ...), plus
  `type_of_algebra()`, which names the kind of single-operator algebra the table supports. A
  table can be built from a table of element indices or a table of element name strings via
  `make_cayley_table`.
- **`FiniteOperator`** (`finite_operator.py`) is a callable wrapper around a `CayleyTable` that
  becomes an algebra's `.op` (or, for a Ring/Field's second table, its multiplication). Called
  with zero arguments it returns the identity (or `None`); with one argument it validates and
  echoes that element; with two it returns their product; with more than two it left-associates,
  e.g. `op(a, b, c, d) == ((a*b)*c)*d`.
- The operator is generically called **"multiplication"**, but for `Group` and its
  single-table ancestors it is more often referred to as the algebra's operation (`op`). Once a
  second operator is introduced (`Ring`/`Field`), the first (inherited) table's operator is
  called **"addition"** (`add`, with additive identity `zero`) and the second table's operator is
  called **"multiplication"** (`mult`, with multiplicative identity `one`, if it exists).
- **`Element`** (`element.py`) wraps a single element string together with its algebra so that
  Python's infix operators (`+`, `-`, `*`, `/`, `**`, unary `-`, and `|` for conjugation) can be
  used for algebraic arithmetic on individual elements. `Element` instances are produced via an
  algebra's `element_map()` method, most conveniently through the `InfixNotation` context
  manager (see [Arithmetic on Elements](#arithmetic-on-elements)).

```python
from finalg.cayley_table import CayleyTable
from finalg.finite_operator import FiniteOperator

table = CayleyTable([[0, 1, 2, 3],
                      [1, 0, 3, 2],
                      [2, 3, 0, 1],
                      [3, 2, 1, 0]])
print(table.is_associative())   # ==> True
print(table.type_of_algebra())  # ==> 'Group'

elements = ('e', 'h', 'v', 'r')
op = FiniteOperator(elements, 'e', table)
print(op('h', 'v'))             # ==> 'r'
```

## Creating Algebras

The recommended (and effectively only) way to build any `FiniteAlgebra` or `FiniteCompositeAlgebra`
is the single factory function `make_finite_algebra` (`make_finite_algebra.py`). Individual class
constructors exist internally but are not part of the intended user interface. `make_finite_algebra`
accepts a variable number of arguments and dispatches accordingly:

- **1 argument** — a path to a JSON file, or a `dict`, holding a serialized `Magma`..`Field`
  algebra (as described below, under [Serialization](#serialization)).
- **4 arguments** — `name, description, elements, table`: builds the most specific of
  `Magma`/`Quasigroup`/`Loop`/`Semigroup`/`Monoid`/`Group` that the table supports.
- **5 arguments**, third and fourth being a `Ring`/`Field` and a `Group`, respectively —
  `ring_or_field, group, operator`: builds a `Module` (if scalars are a plain `Ring`) or a
  `VectorSpace` (if scalars are a `Field`).
- **5 arguments**, otherwise — `name, description, elements, table, table2`: builds a `Ring` or
  `Field`, depending on whether the elements (minus the additive identity) form a commutative
  group under the second table's operator.
- **6 arguments** — as above, plus a `conj_map` dictionary; used internally when constructing
  Cayley–Dickson algebras with explicit conjugation.

```python
from finalg import make_finite_algebra

V4 = make_finite_algebra(
    'V4', 'Klein-4 group',
    ['e', 'h', 'v', 'r'],
    [[0, 1, 2, 3],
     [1, 0, 3, 2],
     [2, 3, 0, 1],
     [3, 2, 1, 0]])
print(type(V4).__name__)  # ==> 'Group'
```

## Serialization

Every single-operator algebra (`Magma` through `Loop`, `Group`) and every `Ring`/`Field` has a
serializable, JSON-compatible representation with these components:

- `name` (`str`)
- `description` (`str`)
- `elements` (list/tuple of `str`) — **all elements must be unique strings**
- `table` (list of lists of `int` or `str`) — an NxN Cayley table for the algebra's operator. By
  default, entries are the *indices* of elements in the `elements` list; entries may instead be
  the element strings themselves (and, if so, only strings from `elements` are allowed — no
  mixing of indices and names).
- `table2` (`Ring`/`Field` only) — a second NxN Cayley table, of the same form, for the
  multiplication operator.
- `conj_map` (`Ring`/`Field` only, optional) — a dictionary mapping each element to its conjugate
  element, used for Cayley–Dickson algebras built with conjugation (see below).

An algebra can be serialized with `to_dict()`, `dumps()` (JSON string), or `dump(filename)`
(JSON file), and reconstructed by passing that dictionary, string, or filename back into
`make_finite_algebra`. `Module` and `VectorSpace` (and their N-dimensional subclasses) do not
currently have a serialized form.

`Module`/`VectorSpace` objects instead have an in-memory representation of five components:
`name`, `description`, a `Ring`/`Field` of scalars, a `Group` of vectors, and a scalar-vector
operator function (see [Modules and Vector Spaces](#modules-and-vector-spaces)).

```python
d = V4.to_dict()
# {'name': 'V4', 'description': 'Klein-4 group', 'elements': ('e', 'h', 'v', 'r'),
#  'table': [[0, 1, 2, 3], [1, 0, 3, 2], [2, 3, 0, 1], [3, 2, 1, 0]]}

json_str = V4.dumps()
V4.dump('v4.json')

V4_copy = make_finite_algebra(d)          # from a dict
V4_copy = make_finite_algebra('v4.json')  # from a JSON file
```

## Arithmetic on Elements

Algebra objects expose their operator as a callable, `.op` (Rings/Fields also expose `.add` and
`.mult`), which can be called postfix-style with two or more elements:

```python
V4.op('h', 'v', V4.inv('r'))   # ==> 'e'
```

For infix arithmetic (`+`, `-`, `*`, `/`, `**`, unary `-`, `|` for conjugation), use the
`InfixNotation` context manager, which yields a dictionary mapping each element name to an
`Element` instance:

```python
from finalg import InfixNotation

with InfixNotation(V4) as v:
    x = v['h'] + v['v'] - v['r']   # ==> Element('e')
```

Which of the infix operators are available depends on what the underlying algebra supports
(e.g., `-` and `/` require `sub`/`div` methods, which only `Group`-and-above or `Field`
algebras provide).

## Algebra-Level Operators: Direct Products and Quotients

Python's `*` and `/` operators are **not** used for arithmetic on individual elements (that's
what `InfixNotation`/`Element` is for). Instead, they're used at the level of whole algebras:

- **`*` (`__mul__`)** — the direct product of two algebras of the same kind (e.g., two `Magma`s,
  or two `Ring`s — mixed kinds raise an error for `Ring.__mul__`). Element names in the product
  are the cross product of the two input element names, joined by a delimiter (`:` by default,
  configurable via `direct_product_delimiter()`), e.g. `V4 * Z2` produces elements like `'e:0'`.
- **`**` (`__pow__`)** — repeated direct product of an algebra with itself, e.g. `V4 ** 3`.
- **`/` (`__truediv__`, `Group` and below in the Ring/Field line)** — the quotient group of a
  `Group` by one of its normal subgroups. Elements of the quotient are representative elements
  of each left coset, prefixed with `~`.

```python
from finalg import generate_cyclic_group

Z2 = generate_cyclic_group(2)
prod = V4 * Z2
print(prod.name, prod.elements)
# ==> V4_x_Z2 ('e:0', 'e:1', 'h:0', 'h:1', 'v:0', 'v:1', 'r:0', 'r:1')

sub = V4.subalgebra_from_elements(['e', 'h'], name='V4_sub')
print(V4.is_normal(sub))          # ==> True (V4 is abelian, so every subgroup is normal)

q = V4 / sub
print(q.name, q.elements)         # ==> V4/V4_sub ('~e', '~v')  (representative elements vary)
```

## What Each Level of the Hierarchy Adds

**`Magma`** (`magma.py`) — the base for anything with a single, closed binary operator. Provides
the broadest set of structural tools, available to every algebra below it:
- direct products (`*`, `**`), element powers (`element_to_power`), reordering elements
  (`reorder_elements`)
- isomorphism checking (`isomorphic`, `is_isomorphic_mapping`, `make_element_mappings`)
- closure, subset/subalgebra discovery (`closure`, `closed_subsets_of_elements`,
  `subalgebra_from_elements`, `proper_subalgebras`)
- generators (`generators`, `is_cyclic`, `generates`, `get_single_generator_set`)
- center (`center`, `center_algebra`)
- left/right cosets (`left_cosets`, `right_cosets`)
- Cayley graphs and diagrams (`make_cayley_graph`, `draw_cayley_diagram`, via NetworkX/Matplotlib)
- a printable summary (`about`)

```python
from finalg import generate_symmetric_group

s3 = generate_symmetric_group(3)
print(s3.is_commutative())        # ==> False
print(s3.center())                # ==> ['(0)(2)']  (just the identity permutation, in this notation)
print(s3.is_cyclic())             # ==> False (needs more than one generator)
print(len(s3.subgroups()))        # ==> 6
```

**`Semigroup`** (`semigroup.py`) — adds associativity (checked at construction), plus
regularity checks (`is_regular`, `weak_inverses`).

**`Monoid`** (`monoid.py`) — adds a required identity element, enabling:
- element order (`element_order`)
- units — elements with a two-sided multiplicative inverse — and the subgroup they form
  (`units`, `units_subgroup`)
- the regular representation of the algebra as (dense or SciPy-sparse) permutation matrices
  (`regular_representation`, `verify_regular_representation`)

```python
mapping, _, elem_to_arr, arr_to_elem = V4.regular_representation()
print(mapping['h'])
# ==> array([[0., 1., 0., 0.],
#            [1., 0., 0., 0.],
#            [0., 0., 0., 1.],
#            [0., 0., 1., 0.]])
```

**`Group`** (`group.py`) — adds required inverses for every element, enabling:
- group subtraction, conjugation, commutators, and the commutator subalgebra (`sub`,
  `conjugate`, `commutator`, `commutators`, `commutator_subalgebra`)
- solvability (`is_solvable`)
- normal subgroup testing, subgroup enumeration (up to isomorphism), and quotient groups
  (`is_normal`, `subgroups`, `unique_proper_subgroups`, `quotient_group`, `/`)

```python
print(V4.is_solvable())           # ==> True (every group of order 4 is solvable)
print(V4.commutators())           # ==> {'e'}  (V4 is abelian, so all commutators are trivial)
```

**`Quasigroup`/`Loop`** (`quasigroup_and_loop.py`) — a separate branch off `Magma` for algebras
with cancellation (unique left/right division) but not necessarily associativity; `Loop` further
requires an identity element. These are thin subclasses today, inheriting all of `Magma`'s
structural methods.

## Rings and Fields

**`Ring`** (`ring.py`) is a `Group` (commutative, by requirement) with a **second** Cayley table
defining an associative multiplication that **distributes** over the inherited addition. The
inherited operator becomes addition (`add`, `zero`/`add_identity`); the new operator is
multiplication (`mult`, `one`/`mult_identity`, if it exists). Notable `Ring` functionality:

- `extract_additive_algebra()` / `extract_multiplicative_algebra()` — pull out the additive
  `Group` or the multiplicative `Semigroup` alone
- `zero_divisors()`, `zero_divisor_pairs()`, `units()` — structural properties of multiplication
- `square_root_mapping()`, `square_roots(elem)` — square roots under multiplication
- `commutator(a, b)` — the ring commutator `(a*b) - (b*a)`
- `element_pairs_where_sum_equals(elem)` / `element_pairs_where_product_equals(elem)`
- direct products of Rings (`*`, requiring both operands be `Ring`s)
- the Cayley–Dickson construction (`make_cayley_dickson_algebra`, `sqr`; see below)

**`Field`** (`field.py`) is a `Ring` whose elements, minus the additive identity, form a
commutative group under multiplication (`is_field` checks this and, if valid, returns that
group as the Field's `mult_abelian_subgroup()`). A trivial one-element Field is not allowed.
Adds `mult_inv(elem)` and `div(x, y)` (`x / y` at the element level, via `Element.__truediv__`).

```python
from finalg import generate_algebra_mod_n

F5 = generate_algebra_mod_n(5)   # 5 is prime ==> Field
print(type(F5).__name__, F5.elements)  # ==> Field ('0', '1', '2', '3', '4')
print(F5.add('2', '3'), F5.mult('2', '3'))  # ==> 0 1
print(F5.units())                # ==> ['1', '2', '3', '4']  (every non-zero element)

R6 = generate_algebra_mod_n(6)   # 6 is not prime ==> Ring
print(type(R6).__name__, R6.zero_divisors())  # ==> Ring ['2', '3', '4']
```

## The Cayley–Dickson Construction

`Ring.make_cayley_dickson_algebra(mu=None, version=1)` builds a new `Ring`/`Field` of twice the
dimension of the original, using elements that are pairs `(a, b)` (named `"a:b"`) and one of
four multiplication conventions, selectable via `version`:

| version | Source | Notes |
|---|---|---|
| 1 (default) | — | `(a,b)(c,d) = (ac − bd, ad + bc)`; no `mu`, no conjugation |
| 2 | Schafer, 1966 | `(a,b)(c,d) = (ac + μ·d·b*, a*·d + c·b)`; requires `mu` |
| 3 | Schafer, 1954 | `(a,b)(c,d) = (ac + μ·d*·b, d·a + b·c*)`; requires `mu` |
| 4 | Baez, 2001 | `(a,b)(c,d) = (ac − d·b*, a*·d + c·b)`; no `mu` |

If `mu` is not supplied, it defaults to the additive inverse of the ring's multiplicative
identity ("−1"), when one exists; versions 2 and 3 require a valid `mu`. The resulting algebra
carries a `conj_map` (accessed via `conjugates()`/`conj(elem)`) so that conjugation and the norm
(`norm(elem) = elem * conj(elem)`) are defined on it. Applying the construction repeatedly to,
e.g., a small finite field yields finite analogues of the complex numbers, quaternions,
octonions, and so on. `sqr()` is an earlier, fixed version of the same idea (equivalent to
version 1). `is_gaussian_prime(elem)` is a helper specific to Cayley–Dickson algebras built over
`generate_algebra_mod_n`.

```python
from finalg import generate_algebra_mod_n

F3 = generate_algebra_mod_n(3)
cda = F3.make_cayley_dickson_algebra()   # version=1, a finite analogue of the complex numbers
print(type(cda).__name__, cda.order)     # ==> Field 9
print(cda.conj('1:2'), cda.norm('1:2'))  # ==> 1:1  2:0
```

## Modules and Vector Spaces

`Module`/`VectorSpace` (`module.py`, `vector_space.py`) represent a different, "composite" kind
of algebra: one with **two** underlying single-operator algebras (a ring/field of scalars and a
group of vectors) tied together by a scalar-vector multiplication function, rather than one or
two Cayley tables of their own. They share the abstract base `FiniteCompositeAlgebra`
(`finite_composite_algebra_ABC.py`), independent of `FiniteAlgebra`. They have no serialized
(JSON) form at present.

A `Module` is built from:
- `ring` — a `Ring` (or `Field`), representing scalars
- `group` — an (abelian) `Group`, representing vectors
- `operator` — a plain Python function/method `(scalar, vector) -> vector`, **not** a
  `FiniteOperator`/`CayleyTable`-backed object

`check_module_conditions` verifies the four module axioms (scaling by one; distributivity of
scalars over vector addition; distributivity of vectors over scalar addition; the associativity
condition relating scalar addition and scaling) hold for the given `ring`, `group`, and
`operator`, and each condition can also be checked individually (`check_scaling_by_one`, etc.).

A **`VectorSpace`** is simply a `Module` whose `ring` is required to be a `Field`.

Two convenience subclasses build an *n*-dimensional module/vector space directly from a single
ring or field, without the caller having to construct the vector group or scalar-vector
operator by hand:

- **`NDimensionalModule(ring, n)`** takes a `Ring`, computes its *n*-fold direct product with
  itself (`ring ** n`) to use as the vector group, and supplies a componentwise scalar-vector
  multiplication (`module_sv_mult`) that splits a vector element like `"a:b:c"` on the direct
  product delimiter, scales each component, and rejoins them.
- **`NDimensionalVectorSpace(field, n)`** does the same, but starting from a `Field`, guaranteeing
  the result is a genuine `VectorSpace`.

Both expose `dimensions`, `origin` (the vector group's identity), and `dot_product(u, v)`
(componentwise scalar multiplication of the two vectors' parts, summed via the scalar ring's
addition). `Module.vector_add(v1, v2)` adds two vectors using the underlying group's operator,
and `Module.about(...)` prints a combined summary of the scalar ring/field and the vector group.

> Note: the *n*-fold direct product of a Field with itself used to build the vectors of an
> `NDimensionalVectorSpace` may itself come out as a `Ring` rather than strictly a `Group` — that's
> fine, since `Ring isa Group`.

```python
from finalg import generate_algebra_mod_n
from finalg.vector_space import NDimensionalVectorSpace

F3 = generate_algebra_mod_n(3)
V = NDimensionalVectorSpace(F3, 2)  # 2-dimensional vector space over F3
print(V.dimensions, V.origin)       # ==> 2 0:0
print(V.vector_add('1:2', '2:1'))   # ==> 0:0
print(V.dot_product('1:2', '2:1'))  # ==> 1
```

## Automatically Generating Algebras

`algebra_generators.py` provides functions (all prefixed `generate_`) for building common
algebras of a given order without hand-writing a Cayley table:

| Function                                                         | Produces                                                                                                                                                                                                                                                                                                                                                    |
|------------------------------------------------------------------|-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| `generate_cyclic_group(order, ...)`                              | Cyclic group `Z_n`                                                                                                                                                                                                                                                                                                                                          |
| `generate_symmetric_group(n, alternating=False, ...)`            | Symmetric group `S_n` (or alternating `A_n`), from permutations                                                                                                                                                                                                                                                                                             |
| `generate_alternating_group(n, ...)`                             | Alternating group `A_n`, from permutations, convenience function                                                                                                                                                                                                                                                                                            |
| `generate_powerset_group(n, ...)`                                | Group on the powerset of `{0,...,n-1}` under symmetric difference                                                                                                                                                                                                                                                                                           |
| `generate_commutative_monoid(order, ...)`                        | Commutative monoid under multiplication mod `n`                                                                                                                                                                                                                                                                                                             |
| `generate_relative_primes_group(n, ...)`                         | Group of integers relatively prime to `n`, under multiplication mod `n`                                                                                                                                                                                                                                                                                     |
| `generate_powerset_ring(n, ...)`                                 | Ring on the powerset of `{0,...,n-1}`: symmetric difference (add) & intersection (mult)                                                                                                                                                                                                                                                                     |
| `generate_algebra_mod_n(n, ...)`                                 | Ring (or, if `n` is prime, Field) of integers mod `n`                                                                                                                                                                                                                                                                                                       |
| `generate_nxn_matrix_algebra(ring, ...)`                         | Ring of 2x2 abstract matrices over a given `Ring` (currently hardcoded to 2x2, though `AbstractMatrix` itself supports arbitrary square shapes)                                                                                                                                                                                                             |
| `generate_dihedral_group(n, ...)`                                | Dihedral group of order `2n`, built from rotation/flip permutation generators                                                                                                                                                                                                                                                                               |
| `generate_algebra_from_element_dict(gen_elem_dict, bin_op, ...)` | General-purpose: closes a dictionary of named generator elements (of *any* type with a defined `bin_op`) under that operation, then builds the corresponding algebra. More flexible than the others, at the cost of requiring more input (a binary operation, an equality test, and functions for making values hashable and for naming generated elements) |

`generate_dihedral_group` is itself implemented on top of `generate_algebra_from_element_dict`,
using `Perm` rotation/flip generators, then relabels elements with `compress_runs` (e.g.
`"ffrrrf"` → `"f^2r^3f"`) for readability.

```python
from finalg import generate_cyclic_group, generate_symmetric_group, generate_dihedral_group

Z6 = generate_cyclic_group(6)
S3 = generate_symmetric_group(3)
D6, elem_dict, iterations = generate_dihedral_group(3)   # dihedral group of order 6
print(Z6.order, S3.order, D6.order)   # ==> 6 6 6
```

## Supporting Modules

- **`permutation.py` — `Perm`.** A standalone permutation implementation (used in place of
  `sympy.combinatorics.permutations`) supporting both 0-based and 1-based permutations,
  composition (`*`, right-to-left), inversion, application to sequences (`Perm(...)`("ABCDE")`),
  parity/sign, and conversion to/from cycle notation (`to_cycles`, `Perm.from_cycles`). See
  `Perm.__init__`'s docstring for the exact cycle-notation conventions, including how singleton
  cycles at the permutation's minimum/maximum value are represented.

  ```python
  from finalg import Perm

  p = Perm((4, 2, 1, 5, 3))
  print(p('ABCDE'))          # ==> 'DBAEC'
  print(p.to_cycles())       # ==> [(1, 4, 5, 3)]
  print(p.inverse().values)  # ==> (3, 2, 5, 1, 4)
  ```

- **`abstract_matrix.py` — `AbstractMatrix`.** Matrices whose entries are elements of a `Ring`
  (or `Field`), supporting `+`, `-`, `*`, scalar multiplication, transpose, minors, determinant
  (via Laplace expansion), cofactor matrix, and inverse (via the adjugate), plus convenience
  constructors `zeros`, `identity`, and `random`.

  ```python
  from finalg import generate_algebra_mod_n
  from finalg.abstract_matrix import AbstractMatrix

  F5 = generate_algebra_mod_n(5)
  m = AbstractMatrix([['1', '2'], ['3', '4']], F5)
  print(m.determinant())  # ==> 3
  ```

- **`about_subalgebras.py`.** Free functions (not methods) for summarizing an algebra's
  subalgebras: `partition_into_isomorphic_lists` groups a list of subalgebras (e.g., subgroups)
  into sublists of mutually isomorphic algebras; `about_subalgebras(alg)` computes an algebra's
  proper subalgebras, partitions them by isomorphism, and prints a readable summary;
  `find_isomorphic_subalgebra` looks up whether a given algebra matches one of the partitions
  found for another algebra.

  ```python
  from finalg import about_subalgebras

  partitions = about_subalgebras(V4)
  # Subalgebras of V4 : Klein-4 group
  #   There is 1 unique proper subalgebra, up to isomorphism, out of 3 total subalgebras.
  #   3 Isomorphic Commutative Normal Groups of order 2 with identity 'e': ...
  ```

- **`examples.py` — `Examples`.** Loads the built-in example algebras shipped as JSON under
  `src/finalg/data/algebras/` (indexed by `examples.json`), accessible as `finalg.examples[i]`
  once the package is imported, with `finalg.examples.about()` listing what's available (dihedral
  and symmetric/alternating groups, the Klein four-group, a quaternion group, small fields and
  rings, a semidihedral group, a quasigroup and loop example, and more).

  ```python
  import finalg

  print(len(finalg.examples))         # ==> 20
  print(finalg.examples[0].name)      # ==> 'A4'
  finalg.examples.about()             # prints the full indexed list
  ```

- **`my_math.py`.** Small number-theoretic helpers not tied to any particular algebra:
  `is_relatively_prime`, `relative_primes`, `totient`, `divisors`, and `xgcd` (extended
  Euclidean algorithm).

  ```python
  from finalg.my_math import totient, xgcd

  print(totient(9))     # ==> 6
  print(xgcd(240, 46))  # ==> (2, -9, 47)
  ```

- **`utilities.py`.** General helpers used across the package: `np_arr_to_tuple` (also usable as
  the `make_immutable` argument to `generate_algebra_from_element_dict`), `get_duplicates`,
  `all_strings`, `powerset`, `compress_runs`, `cayley_graph_to_json` (NetworkX node-link
  serialization of a Cayley graph), and `make_table_from_xml` (a helper for converting tables
  copied from the "Groupprops" wiki into `finalg`'s list-of-lists format).

  ```python
  from finalg.utilities import compress_runs

  print(compress_runs("ffrrrf"))  # ==> 'f^2r^3f'
  ```

## Repository Layout

```
src/finalg/
├── finite_algebra_ABC.py           # FiniteAlgebra (ABC)
├── finite_composite_algebra_ABC.py # FiniteCompositeAlgebra (ABC)
├── magma.py                        # Magma
├── semigroup.py                    # Semigroup
├── monoid.py                       # Monoid
├── group.py                        # Group
├── quasigroup_and_loop.py          # Quasigroup, Loop
├── ring.py                         # Ring
├── field.py                        # Field
├── module.py                       # Module, NDimensionalModule
├── vector_space.py                 # VectorSpace, NDimensionalVectorSpace
├── cayley_table.py                 # CayleyTable, make_cayley_table, about_tables
├── finite_operator.py              # FiniteOperator
├── element.py                      # Element
├── infix_notation.py               # InfixNotation
├── make_finite_algebra.py          # make_finite_algebra (the factory function)
├── algebra_generators.py           # generate_* functions
├── about_subalgebras.py            # subalgebra/isomorphism summary tools
├── abstract_matrix.py              # AbstractMatrix
├── permutation.py                  # Perm
├── examples.py                     # Examples (built-in example algebras)
├── my_math.py                      # number-theory helpers
├── utilities.py                    # general helpers
└── data/algebras/                  # JSON-serialized example algebras
tests/                              # pytest test suite
docs/                                # Sphinx documentation (finalg.readthedocs.io)
```

`finalg` requires Python ≥3.10 and depends on NumPy, SciPy, SymPy, NetworkX, and Matplotlib.

## References

- Pinter, C.C., *A Book of Abstract Algebra*, McGraw-Hill, 1982.
- Sawyer, W.W., *A Concrete Approach to Abstract Algebra*, Dover Publications, 1978.
- Carter, N.C., *Visual Group Theory*, Mathematical Association of America, 2009.
- W. B. Vasantha Kandasamy, *Groupoids and Smarandache Groupoids*.
- Schafer, R.D., 1954 and 1966 treatments of the Cayley–Dickson construction.
- Baez, J., "The Octonions", 2001, for the Cayley–Dickson conventions used in `version=4`.

See `docs/90_resources.rst` for a fuller list of books, papers, and websites that informed this
project.
