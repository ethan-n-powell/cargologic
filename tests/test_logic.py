"""Logic layer: parsing and minimisation correctness."""

import pytest

from cargologic.logic import (
    Implicant,
    LogicError,
    all_functions,
    minimise,
    parse_expression,
    parse_table,
    prime_implicants,
)


# ---------------------------------------------------------------------------- parsing


@pytest.mark.parametrize(
    "expr,ones",
    [
        ("A & B", {3}),
        ("A and B", {3}),
        ("A | B", {1, 2, 3}),
        ("A or B", {1, 2, 3}),
        ("A ^ B", {1, 2}),
        ("A xor B", {1, 2}),
        ("!A", {0}),
        ("not A", {0}),
        ("~A", {0}),
        ("A & !B", {2}),
        ("!(A & B)", {0, 1, 2}),
        ("!(A | B)", {0}),
        ("A B", {3}),  # implicit AND
        ("(A | B) & !C", {2, 4, 6}),
    ],
)
def test_parse_expression(expr, ones):
    assert parse_expression(expr).ones == ones


@pytest.mark.parametrize("gate,ones", [("AND", {3}), ("XOR", {1, 2}), ("NAND", {0, 1, 2})])
def test_named_gates(gate, ones):
    assert parse_expression(gate).ones == ones


def test_input_order_is_msb_first():
    """Variable 0 is the most significant bit, so '10' means A=1, B=0."""
    f = parse_expression("A & !B")
    assert f.ones == {0b10}
    assert f.inputs == ["A", "B"]


def test_explicit_inputs_fixes_order():
    f = parse_expression("B", inputs=["A", "B"])
    assert f.n == 2
    assert f.ones == {0b01, 0b11}


def test_parse_expression_rejects_undeclared_variable():
    with pytest.raises(LogicError):
        parse_expression("A & Z", inputs=["A", "B"])


@pytest.mark.parametrize("bad", ["", "A &", "A & & B", "(A", "A $ B"])
def test_parse_expression_errors(bad):
    with pytest.raises(LogicError):
        parse_expression(bad)


@pytest.mark.parametrize(
    "spec,ones",
    [
        ("00:0,01:1,10:1,11:0", {1, 2}),
        ("0110", {1, 2}),
        ("1,2", {1, 2}),
        ("0001", {3}),
    ],
)
def test_parse_table(spec, ones):
    assert parse_table(spec).ones == ones


def test_parse_table_dont_cares():
    f = parse_table("00:1,01:1,10:x,11:0")
    assert f.ones == {0, 1}
    assert f.dont_cares == {2}


def test_parse_table_rejects_decimal_row_labels_with_guidance():
    with pytest.raises(LogicError, match="not binary"):
        parse_table("0:1,1:1,2:x,3:0")


def test_parse_table_rejects_non_power_of_two_column():
    with pytest.raises(LogicError, match="power of two"):
        parse_table("011")


# ----------------------------------------------------------------------- minimisation


def test_minimise_and():
    f = parse_expression("A & B")
    cov = f.minimise()
    assert cov.n_terms == 1
    assert cov.terms[0].positive_literals() == [0, 1]
    assert cov.terms[0].negative_literals() == []


def test_minimise_collapses_or_to_single_literals():
    cov = parse_expression("A | B").minimise()
    assert cov.n_terms == 2
    # Each term should be a single positive literal, not a full minterm.
    assert all(t.n_fixed == 1 for t in cov.terms)


def test_minimise_exploits_dont_cares():
    """A don't-care should let the cover shrink."""
    strict = parse_table("1000").minimise()        # only !A!B
    relaxed = parse_table("10x0").minimise()       # !A!B, with !A B free
    assert relaxed.n_literals <= strict.n_literals


def test_minimise_xor_needs_two_terms():
    cov = parse_expression("A ^ B").minimise()
    assert cov.n_terms == 2


def test_minimise_majority_needs_three_terms():
    cov = parse_expression("MAJORITY3").minimise()
    assert cov.n_terms == 3


def test_minimise_tautology_is_one_ungated_term():
    cov = parse_table("1111").minimise()
    assert cov.n_terms == 1
    assert cov.terms[0].is_tautology()


def test_minimise_contradiction_is_empty():
    assert parse_table("0000").minimise().n_terms == 0


@pytest.mark.parametrize("n", [1, 2, 3])
def test_every_function_cover_reproduces_its_truth_table(n):
    """The heart of it: minimisation must never change the function.

    All 2, 16 and 256 functions of 1, 2 and 3 inputs.
    """
    for f in all_functions(n):
        cov = f.minimise()
        for m in range(1 << n):
            assert cov.evaluate(m) == f.evaluate(m), (f.name, m)


def test_four_input_functions_sample():
    """Spot-check 4 inputs, where exhaustive enumeration would be 65536 functions."""
    import random

    rng = random.Random(0)
    for _ in range(200):
        ones = {m for m in range(16) if rng.random() < 0.5}
        cov = minimise(4, ones)
        for m in range(16):
            assert cov.evaluate(m) == (m in ones)


def test_prime_implicants_ignore_pure_dont_care_cubes():
    primes = prime_implicants(2, ones={3}, dont_cares={0, 1, 2})
    assert all(any(m == 3 for m in p.minterms()) for p in primes)


def test_cover_cost_key_prefers_fewer_terms_then_fewer_blockers():
    """Negative literals are nearly free; positive literals each cost a blocker."""
    nand = parse_expression("!(A & B)").minimise()
    orr = parse_expression("A | B").minimise()
    assert nand.n_terms == orr.n_terms == 2
    # NAND minimises to !B | !A -- all negative, so no blockers at all.
    assert nand.n_positive_literals == 0
    assert orr.n_positive_literals == 2
    assert nand.cost_key < orr.cost_key


# --------------------------------------------------------------------------- implicant


def test_implicant_pattern_and_literals():
    imp = Implicant(bits=0b10, mask=0b11, n=2)
    assert imp.pattern() == "10"
    assert imp.positive_literals() == [0]
    assert imp.negative_literals() == [1]
    assert imp.minterms() == [0b10]


def test_implicant_with_free_variable_covers_both():
    imp = Implicant(bits=0b10, mask=0b10, n=2)  # A=1, B free
    assert sorted(imp.minterms()) == [0b10, 0b11]
    assert imp.to_str(["A", "B"]) == "A"
