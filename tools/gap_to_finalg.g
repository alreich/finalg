# gap_to_finalg.g  --  write a GAP group as a finalg JSON file (pretty-printed).
#
# Usage, in GAP:
#
#     gap> Read("gap_to_finalg.g");
#     gap> FinalgJSON(SmallGroup(32, 8), "SG32_8",
#                     "GAP SmallGroup(32,8): C2 . ((C4 x C2) : C2)", "G32_08_SG32_8.json");
#
# The file can then be read by finalg:  make_finite_algebra("G32_08_SG32_8.json").
#
# Elements of a solvable group are named by their polycyclic exponent vectors with respect to GAP's
# pcgs g1, g2, ..., gr: the identity is "e", and an element g1^a1 g2^a2 ... gr^ar is written with the
# factors that occur, e.g. "g1g3" or "g2^2g4" (separated by "*" if r >= 10). Other groups get the names
# "g0", "g1", ... in the order of GAP's Elements(G).  The identity is always the first element.

FinalgJSON := function(G, name, description, filename)
  local pcgs, vecs, els, names, n, index, table, row, i, j, vec, parts, sep, quote, wrap, out, width, line, piece, lines, rows;

  quote := function(s)
    local r, c;
    r := "\"";
    for c in s do
      if c = '"' or c = '\\' then Add(r, '\\'); fi;
      Add(r, c);
    od;
    Add(r, '"');
    return r;
  end;

  if IsSolvableGroup(G) then
    pcgs := Pcgs(G);
    vecs := Cartesian(List(RelativeOrders(pcgs), r -> [0..r-1]));    # lexicographic; the identity comes first
    els := List(vecs, v -> PcElementByExponents(pcgs, v));
    sep := "";
    if Length(pcgs) >= 10 then sep := "*"; fi;
    names := [];
    for vec in vecs do
      parts := [];
      for i in [1..Length(vec)] do
        if vec[i] = 1 then Add(parts, Concatenation("g", String(i)));
        elif vec[i] > 1 then Add(parts, Concatenation("g", String(i), "^", String(vec[i])));
        fi;
      od;
      if Length(parts) = 0 then Add(names, "e"); else Add(names, JoinStringsWithSeparator(parts, sep)); fi;
    od;
    n := Length(els);
    index := a -> PositionSorted(vecs, ExponentsOfPcElement(pcgs, a));
  else
    els := Concatenation([Identity(G)], Filtered(AsSSortedList(G), g -> g <> Identity(G)));
    names := List([0..Length(els) - 1], i -> Concatenation("g", String(i)));
    n := Length(els);
    index := a -> Position(els, a);
  fi;

  table := List(els, a -> List(els, b -> index(a * b) - 1));

  # Wrap the element names to 100 columns, like finalg's other pretty-printed files
  width := 100;
  lines := [];
  line := "";
  for piece in List(names, quote) do
    if Length(line) > 0 and 4 + Length(line) + 2 + Length(piece) + 1 > width then
      Add(lines, Concatenation(line, ","));
      line := piece;
    elif Length(line) > 0 then
      line := Concatenation(line, ", ", piece);
    else
      line := piece;
    fi;
  od;
  Add(lines, line);

  out := OutputTextFile(filename, false);
  SetPrintFormattingStatus(out, false);
  AppendTo(out, "{\n");
  AppendTo(out, "  \"name\": ", quote(name), ",\n");
  AppendTo(out, "  \"description\": ", quote(description), ",\n");
  AppendTo(out, "  \"elements\": [\n");
  for line in lines do AppendTo(out, "    ", line, "\n"); od;
  AppendTo(out, "  ],\n");
  AppendTo(out, "  \"table\": [\n");
  for i in [1..n] do
    AppendTo(out, "    [", JoinStringsWithSeparator(List(table[i], String), ", "), "]");
    if i < n then AppendTo(out, ",\n"); else AppendTo(out, "\n"); fi;
  od;
  AppendTo(out, "  ]\n}\n");
  CloseStream(out);
  return n;
end;
