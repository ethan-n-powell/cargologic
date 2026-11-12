"""Architecture, delivery, synthesis and end-to-end compilation."""

import pytest

import cargologic as cl
from cargologic.architecture import can_share_payload
from cargologic.delivery import (
    lentiviral_relative_titre,
    minimum_vectors_needed,
    vector_catalogue,
)
from cargologic.logic import all_functions, parse_expression
from cargologic.parts import AAV_ITR_BP, Parts
from cargologic.synthesis import CircuitSpec, PartChoice, SynthesisError, build_design


def test_size_is_sum_of_elements():
    r = cl.compile_circuit("A & !B", payload=1000, optimise=False)
    assert r.design.total_bp == sum(e.bp for e in r.design.elements)


def test_feature_table_is_contiguous_and_covers_the_whole_design():
    r = cl.compile_circuit("A | B", payload=800, optimise=False)
    rows = r.design.feature_table()
    assert rows[0]["start"] == 1
    for a, b in zip(rows, rows[1:]):
        assert b["start"] == a["end"] + 1, "features must tile without gaps or overlap"
    assert rows[-1]["end"] == r.design.total_bp


def test_payload_size_shows_up_one_for_one():
    """Growing the payload by N bp grows a single-cassette design by exactly N bp."""
    a = cl.compile_circuit("A & B", payload=1000, optimise=False)
    b = cl.compile_circuit("A & B", payload=1500, optimise=False)
    assert b.design.total_bp - a.design.total_bp == 500


def test_extra_product_term_costs_more_than_a_part_swap():
    one_term = cl.compile_circuit("A & B", payload=1500, optimise=False)
    two_term = cl.compile_circuit("A ^ B", payload=1500, optimise=False)
    term_cost = two_term.design.total_bp - one_term.design.total_bp

    biggest = cl.compile_circuit("A & B", payload=1500, optimise=False, promoter="CAG")
    smallest = cl.compile_circuit("A & B", payload=1500, optimise=False, promoter="minCMV")
    promoter_span = biggest.design.total_bp - smallest.design.total_bp

    assert term_cost > promoter_span


def test_or_can_share_payload_but_xor_cannot():
    assert can_share_payload(parse_expression("A | B").minimise())
    assert not can_share_payload(parse_expression("A ^ B").minimise())


def test_sharing_makes_or_smaller_than_xor_at_equal_term_count():
    orr = cl.compile_circuit("A | B", payload=1500)
    xor = cl.compile_circuit("A ^ B", payload=1500)
    assert orr.covers["out"].n_terms == xor.covers["out"].n_terms == 2
    assert orr.design.strategy == "tandem_promoters"
    assert xor.design.strategy == "independent"
    assert orr.circuit_bp < xor.circuit_bp


def test_identical_logic_on_two_outputs_shares_one_cassette():
    r = cl.compile_circuit({"x": "A & !B", "y": "A & !B"}, payload={"x": 1000, "y": 1000})
    assert r.design.n_terms == 1
    assert any("P2A" in n for n in r.design.notes)


def test_forcing_tandem_on_unshareable_cover_raises():
    f = parse_expression("A ^ B")
    spec = CircuitSpec(
        inputs=f.inputs, functions=[f], payloads={"out": Parts().resolve_payload(500)}
    )
    with pytest.raises(SynthesisError, match="negative-literal"):
        build_design(spec, Parts(), PartChoice(), force_strategy="tandem_promoters")


def test_aav_single_budget_is_capacity_minus_itrs():
    assert vector_catalogue()["aav_single"].budget() == 4700 - 2 * AAV_ITR_BP


def test_oversized_aav_is_flagged_over_not_merely_tight():
    verdict = vector_catalogue()["aav_single"].evaluate(6000)
    assert not verdict.fits
    assert verdict.severity == "over"
    assert verdict.headroom_bp < 0


def test_near_limit_is_flagged_tight():
    v = vector_catalogue()["aav_single"]
    verdict = v.evaluate(int(0.97 * v.budget()))
    assert verdict.fits
    assert verdict.severity == "tight"


def test_mrna_cannot_carry_a_dna_circuit():
    verdict = vector_catalogue()["mrna_lnp"].evaluate(1000)
    assert not verdict.fits
    assert verdict.severity == "infeasible"


def test_lentiviral_titre_is_flat_below_reference_then_falls():
    assert lentiviral_relative_titre(2000) == pytest.approx(1.0)
    assert lentiviral_relative_titre(4400) == pytest.approx(1.0)
    assert lentiviral_relative_titre(8000) < 0.5
    # Anchored on ~2 orders of magnitude lost between 4.4 kb and 11.1 kb.
    assert lentiviral_relative_titre(11100) == pytest.approx(0.01, rel=0.05)


def test_titre_is_monotonically_decreasing():
    titres = [lentiviral_relative_titre(s) for s in range(4000, 14000, 500)]
    assert all(a >= b for a, b in zip(titres, titres[1:]))


@pytest.mark.parametrize(
    "bp,expected", [(1000, 1), (4440, 1), (4441, 2), (8000, 2), (9000, 3)]
)
def test_minimum_vector_count(bp, expected):
    assert minimum_vectors_needed(bp) == expected


def test_single_vector_count_agrees_with_fit_verdict():
    v = vector_catalogue()["aav_single"]
    for bp in (1000, 3000, 4400, 4440):
        assert (minimum_vectors_needed(bp) == 1) == v.evaluate(bp).fits


def test_search_does_not_waste_budget_for_unresolvable_gains():
    """Given equivalent dynamic range, the compiler must take the smaller design."""
    r = cl.compile_circuit("A & !B", payload="CAR_2g", n_samples=300)
    front = r.pareto
    assert front
    assert r.circuit_bp <= min(c.circuit_bp for c in front) + 1


def test_pareto_front_is_not_degenerate():
    """Designs with model-indistinguishable dynamic range collapse to one point."""
    r = cl.compile_circuit("A & !B", payload="CAR_2g", n_samples=300)
    drs = [c.dynamic_range for c in r.pareto]
    for a, b in zip(drs, drs[1:]):
        assert b > a * 1.05, "front should not list equivalent-performing designs"


def test_explicit_parts_are_respected():
    r = cl.compile_circuit("A & B", payload=500, promoter="CMV", terminator="SPA")
    assert r.choice.promoter == "CMV"
    assert r.choice.terminator == "SPA"


def test_not_gate_is_excision_limited_not_terminator_limited():
    rep = cl.compile_circuit("A & !B", payload=1500, n_samples=800).leakage["out"]
    assert rep.limiting is not None
    assert rep.limiting.name == "incomplete excision"
    # "11" is the OFF state with the NOT input high.
    assert rep.worst_off_state == "11"


def test_improving_excision_helps_a_not_gate_but_better_terminators_do_not():
    base = cl.compile_circuit("A & !B", payload=1500, n_samples=800).leakage["out"]

    better_rec = Parts()
    better_rec.override_recombinase("FlpO", excision_efficiency=(0.98, 0.999))
    improved = cl.compile_circuit(
        "A & !B", payload=1500, parts=better_rec, n_samples=800
    ).leakage["out"]
    assert improved.dynamic_range_median > 1.5 * base.dynamic_range_median

    tighter_term = cl.compile_circuit(
        "A & !B", payload=1500, optimise=False, blocker="tandem", terminator="tandem",
        n_samples=800,
    ).leakage["out"]
    loose_term = cl.compile_circuit(
        "A & !B", payload=1500, optimise=False, blocker="SPA", terminator="SPA",
        n_samples=800,
    ).leakage["out"]
    assert tighter_term.dynamic_range_median == pytest.approx(
        loose_term.dynamic_range_median, rel=0.15
    )


def test_pure_and_gate_has_no_off_state_requiring_excision():
    """AND uses only blockers, so its OFF states are terminator-limited."""
    rep = cl.compile_circuit("A & B", payload=1500, n_samples=800).leakage["out"]
    assert rep.limiting is not None
    assert rep.limiting.name == "terminator read-through"


def test_constant_on_has_no_off_state():
    r = cl.compile_circuit("1", payload=500, inputs=["A", "B"], n_samples=100)
    assert r.leakage["out"].dynamic_range_median == float("inf")


def test_constant_off_emits_no_cassette():
    r = cl.compile_circuit("0000", payload=500, as_table=True, n_samples=100)
    assert r.design.n_terms == 0
    assert r.design.total_bp == 0
    assert "never expressed" in " ".join(r.leakage["out"].warnings)


def test_stacked_blockers_warn_about_weak_on():
    r = cl.compile_circuit("A & B & C", payload=800, n_samples=400)
    assert any("blockers in series" in w for w in r.leakage["out"].warnings)


def test_predictions_are_deterministic_for_a_fixed_seed():
    a = cl.compile_circuit("A & !B", payload=1200, n_samples=400)
    b = cl.compile_circuit("A & !B", payload=1200, n_samples=400)
    assert a.leakage["out"].dynamic_range_median == b.leakage["out"].dynamic_range_median


def test_all_sixteen_two_input_functions_compile():
    for f in all_functions(2):
        column = "".join("1" if m in f.ones else "0" for m in range(4))
        r = cl.compile_circuit(
            column, payload="CAR_2g", as_table=True, inputs=f.inputs, n_samples=100
        )
        assert r.feasible, f.name
        for m in range(4):
            assert r.covers["out"].evaluate(m) == f.evaluate(m)


def test_cas9_does_not_fit_one_aav_but_fits_two():
    r = cl.compile_circuit("A & B", payload="SpCas9", vector="aav_single", n_samples=100)
    assert not r.feasible
    assert cl.minimum_aav_count(r) == 2
    dual = next(v for v in r.all_verdicts if "dual" in v.vector)
    assert dual.fits


def test_carrying_own_recombinases_is_accounted_for():
    plain = cl.compile_circuit("A & !B", payload=1500, n_samples=100)
    loaded = cl.compile_circuit(
        "A & !B", payload=1500, include_recombinases=True, n_samples=100
    )
    assert loaded.design.recombinase_bp > 2000
    assert loaded.circuit_bp > plain.circuit_bp


def test_three_inputs_get_three_orthogonal_recombinases():
    r = cl.compile_circuit("A & B & !C", payload=1000, n_samples=100)
    names = [x.name for x in r.design.recombinases]
    assert len(names) == len(set(names)) == 3


def test_too_many_inputs_raises_rather_than_inventing_recombinases():
    with pytest.raises(ValueError, match="orthogonal"):
        cl.compile_circuit("A & B & C & D & E & F & G", payload=500, n_samples=50)


def test_crosstalking_recombinases_are_not_auto_selected_together():
    chosen = Parts().orthogonal_set(5)
    assert not ("Bxb1" in chosen and "phiC31" in chosen)


def test_report_renders_and_names_the_failing_state():
    text = cl.compile_circuit("A & !B", payload="CAR_2g", n_samples=400).report()
    assert "FITS" in text
    assert "LOGIC MINIMISATION" in text
    assert "what is limiting" in text
    assert "CAVEATS" in text


def test_genbank_is_parseable_and_length_consistent():
    r = cl.compile_circuit("A | B", payload=900, n_samples=100)
    gb = r.to_genbank("test_design")
    assert gb.startswith("LOCUS")
    assert gb.rstrip().endswith("//")
    assert f"{r.design.total_bp} bp" in gb
    seq = "".join(
        part
        for line in gb.split("ORIGIN")[1].splitlines()
        for part in line.split()
        if set(part) <= set("n")
    )
    assert len(seq) == r.design.total_bp


def test_unknown_vector_and_payload_errors_are_clear():
    with pytest.raises(KeyError, match="unknown vector"):
        cl.compile_circuit("A & B", vector="banana")
    with pytest.raises(KeyError, match="unknown payload"):
        cl.compile_circuit("A & B", payload="not_a_real_payload")


def test_multi_output_missing_payload_raises():
    with pytest.raises(SynthesisError, match="payload"):
        cl.compile_circuit({"x": "A", "y": "B"}, payload={"x": 500})


def test_parts_override_round_trips():
    p = Parts()
    p.override_recombinase("Cre", excision_efficiency=(0.5, 0.6))
    assert p.recombinases["Cre"].excision_efficiency == (0.5, 0.6)
    assert p.recombinases["Cre"].site_bp == 34  # untouched


def test_not_gates_cannot_be_rescued_by_better_terminators_but_and_gates_can():
    """An AND OFF state is held by a terminator, a NOT OFF state by an excision event."""
    and_gate = cl.compile_circuit("A & B", payload="CAR_2g", n_samples=600)
    not_gate = cl.compile_circuit("A & B & !C", payload="CAR_2g", n_samples=600)

    assert and_gate.leakage["out"].limiting.name == "terminator read-through"
    assert not_gate.leakage["out"].limiting.name == "incomplete excision"

    assert len(and_gate.pareto) >= 3
    best, worst = max(and_gate.pareto, key=lambda c: c.dynamic_range), min(
        and_gate.pareto, key=lambda c: c.dynamic_range
    )
    assert best.dynamic_range > 2 * worst.dynamic_range
    assert best.circuit_bp > worst.circuit_bp

    assert len(not_gate.pareto) == 1
