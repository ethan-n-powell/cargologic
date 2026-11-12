"""Payload-constrained logic compilation for mammalian recombinase circuits."""

from .compiler import CompileResult, compile_circuit, minimum_aav_count
from .delivery import (
    DeliveryVerdict,
    evaluate_all,
    lentiviral_relative_titre,
    vector_catalogue,
)
from .logic import (
    BooleanFunction,
    Cover,
    Implicant,
    LogicError,
    all_functions,
    minimise,
    parse_expression,
    parse_table,
)
from .parts import Parts, Payload
from .synthesis import (
    Candidate,
    CircuitSpec,
    PartChoice,
    SynthesisError,
    build_design,
    pareto_front,
    search,
)

__version__ = "0.1.0"

__all__ = [
    "compile_circuit",
    "CompileResult",
    "minimum_aav_count",
    "Parts",
    "Payload",
    "parse_expression",
    "parse_table",
    "all_functions",
    "minimise",
    "BooleanFunction",
    "Cover",
    "Implicant",
    "LogicError",
    "SynthesisError",
    "CircuitSpec",
    "PartChoice",
    "Candidate",
    "build_design",
    "search",
    "pareto_front",
    "vector_catalogue",
    "evaluate_all",
    "lentiviral_relative_titre",
    "DeliveryVerdict",
    "__version__",
]
