"""cargologic -- payload-constrained logic compilation for mammalian recombinase circuits.

You give it a truth table and a delivery format. It gives you a DNA architecture, an
exact size, a packaging verdict, and a modelled dynamic range with uncertainty.

    >>> import cargologic as cl
    >>> r = cl.compile_circuit("A and not B", payload="CAR_2g", vector="aav_single")
    >>> r.feasible
    True

The size verdict is arithmetic over physical packaging limits and is trustworthy. The
expression predictions are a model and are not. The report keeps the two apart on
purpose.
"""

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
