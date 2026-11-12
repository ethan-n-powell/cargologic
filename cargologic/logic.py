"""Boolean specification and two-level logic minimisation."""

from __future__ import annotations

import re
from dataclasses import dataclass


class LogicError(ValueError):
    """Raised for malformed Boolean specifications."""


@dataclass(frozen=True, order=True)
class Implicant:
    """A product term over n variables. mask marks the constrained variables and bits
    their values; variable i is bit (n - 1 - i), so variable 0 is most significant.
    """

    bits: int
    mask: int
    n: int

    def covers(self, minterm: int) -> bool:
        return (minterm & self.mask) == self.bits

    @property
    def n_fixed(self) -> int:
        return bin(self.mask).count("1")

    def _bit_of(self, i: int) -> int:
        return 1 << (self.n - 1 - i)

    def literal(self, i: int) -> int | None:
        """1 if variable i must be true, 0 if false, None if it is a don't-care."""
        b = self._bit_of(i)
        if not (self.mask & b):
            return None
        return 1 if (self.bits & b) else 0

    def positive_literals(self) -> list[int]:
        return [i for i in range(self.n) if self.literal(i) == 1]

    def negative_literals(self) -> list[int]:
        return [i for i in range(self.n) if self.literal(i) == 0]

    def minterms(self) -> list[int]:
        out = [self.bits]
        for b in (self._bit_of(i) for i in range(self.n) if self.literal(i) is None):
            out = out + [v | b for v in out]
        return sorted(out)

    def is_tautology(self) -> bool:
        return self.mask == 0

    def to_str(self, names: list[str]) -> str:
        parts = []
        for i, nm in enumerate(names):
            lit = self.literal(i)
            if lit is not None:
                parts.append(nm if lit else f"!{nm}")
        return " & ".join(parts) if parts else "1"

    def pattern(self, names: list[str] | None = None) -> str:
        """Compact '1-0' style pattern."""
        return "".join(
            "-" if (lit := self.literal(i)) is None else str(lit) for i in range(self.n)
        )


def _combine(a: Implicant, b: Implicant) -> Implicant | None:
    """Merge two cubes that differ in exactly one fixed variable."""
    if a.mask != b.mask:
        return None
    diff = a.bits ^ b.bits
    if diff == 0 or (diff & (diff - 1)) != 0:
        return None
    new_mask = a.mask & ~diff
    return Implicant(a.bits & new_mask, new_mask, a.n)


def prime_implicants(
    n: int, ones: set[int], dont_cares: set[int] = frozenset()
) -> list[Implicant]:
    """Quine-McCluskey prime implicant generation."""
    terms = set(ones) | set(dont_cares)
    if not terms:
        return []
    full_mask = (1 << n) - 1
    current = {Implicant(t, full_mask, n) for t in terms}
    primes: set[Implicant] = set()

    while current:
        used: set[Implicant] = set()
        nxt: set[Implicant] = set()
        # Grouped by set-bit count, so only adjacent groups need comparing.
        by_weight: dict[int, list[Implicant]] = {}
        for imp in sorted(current):
            by_weight.setdefault(bin(imp.bits & imp.mask).count("1"), []).append(imp)
        for w in sorted(by_weight):
            for a in by_weight[w]:
                for b in by_weight.get(w + 1, []):
                    merged = _combine(a, b)
                    if merged is not None:
                        nxt.add(merged)
                        used.update((a, b))
        primes |= current - used
        current = nxt

    return sorted(p for p in primes if any(m in ones for m in p.minterms()))


@dataclass
class Cover:
    """A chosen set of product terms covering a function's ON-set."""

    terms: list[Implicant]

    @property
    def n_terms(self) -> int:
        return len(self.terms)

    @property
    def n_positive_literals(self) -> int:
        return sum(len(t.positive_literals()) for t in self.terms)

    @property
    def n_literals(self) -> int:
        return sum(t.n_fixed for t in self.terms)

    @property
    def cost_key(self) -> tuple[int, int, int]:
        """Sort key: fewer cassettes, then fewer blockers, then fewer sites."""
        return (self.n_terms, self.n_positive_literals, self.n_literals)

    def evaluate(self, minterm: int) -> bool:
        return any(t.covers(minterm) for t in self.terms)

    def to_str(self, names: list[str]) -> str:
        return "  |  ".join(t.to_str(names) for t in self.terms) if self.terms else "0"


def minimise(
    n: int,
    ones: set[int],
    dont_cares: set[int] = frozenset(),
    exact_limit: int = 10,
) -> Cover:
    """Minimal sum-of-products cover; greedy only above exact_limit."""
    ones = set(ones)
    if not ones:
        return Cover([])
    primes = prime_implicants(n, ones, dont_cares)
    if not primes:
        return Cover([])

    for p in primes:
        if p.is_tautology():
            return Cover([p])

    coverage: dict[int, list[int]] = {}
    for m in sorted(ones):
        hit = [i for i, p in enumerate(primes) if p.covers(m)]
        if not hit:  # pragma: no cover
            raise LogicError(f"minterm {m} uncovered by any prime implicant")
        coverage[m] = hit

    # Essential primes: the sole cover of some minterm.
    essential = {hits[0] for hits in coverage.values() if len(hits) == 1}
    covered = {m for i in essential for m in primes[i].minterms()}
    remaining = ones - covered

    if not remaining:
        return Cover(sorted((primes[i] for i in essential), key=lambda p: p.pattern()))

    sub_coverage = {m: coverage[m] for m in sorted(remaining)}
    if max(len(v) for v in sub_coverage.values()) <= exact_limit:
        choices = _petrick(sub_coverage)
    else:  # pragma: no cover
        choices = [_greedy(sub_coverage, primes)]

    candidates = [
        Cover(
            sorted(
                (primes[i] for i in sorted(essential | extra)),
                key=lambda p: p.pattern(),
            )
        )
        for extra in choices
    ]
    return min(candidates, key=lambda c: c.cost_key)


def _petrick(coverage: dict[int, list[int]]) -> list[set[int]]:
    """Expand the product-of-sums cover condition, absorbing supersets as it goes."""
    solutions: list[frozenset[int]] = [frozenset()]
    for hits in coverage.values():
        grown: list[frozenset[int]] = []
        for sol in solutions:
            if sol & set(hits):
                grown.append(sol)
            else:
                grown += [sol | {h} for h in hits]
        minimal: list[frozenset[int]] = []
        for s in sorted(set(grown), key=len):
            if not any(m <= s for m in minimal):
                minimal.append(s)
        solutions = minimal
    smallest = min(len(s) for s in solutions)
    return [set(s) for s in solutions if len(s) == smallest]


def _greedy(
    coverage: dict[int, list[int]], primes: list[Implicant]
) -> set[int]:  # pragma: no cover
    """Greedy set cover fallback for pathologically large functions."""
    need = set(coverage)
    chosen: set[int] = set()
    while need:
        best_i, best_gain = None, -1
        for i, p in enumerate(primes):
            gain = len(need & set(p.minterms()))
            if gain > best_gain:
                best_i, best_gain = i, gain
        if best_i is None or best_gain <= 0:
            break
        chosen.add(best_i)
        need -= set(primes[best_i].minterms())
    return chosen


_TOKEN_RE = re.compile(
    r"\s*(?:(?P<lparen>\()|(?P<rparen>\))|(?P<op>&&|\|\||\^|&|\||~|!)"
    r"|(?P<word>[A-Za-z_][A-Za-z_0-9]*)|(?P<const>[01])|(?P<bad>\S))"
)

_WORD_OPS = {"and": "&", "or": "|", "not": "!", "xor": "^"}


def _tokenise(expr: str) -> list[tuple[str, str]]:
    toks: list[tuple[str, str]] = []
    pos = 0
    while pos < len(expr):
        m = _TOKEN_RE.match(expr, pos)
        if not m or m.end() == pos:
            break
        pos = m.end()
        if m.group("lparen"):
            toks.append(("lparen", "("))
        elif m.group("rparen"):
            toks.append(("rparen", ")"))
        elif op := m.group("op"):
            toks.append(("op", {"&&": "&", "||": "|", "~": "!"}.get(op, op)))
        elif w := m.group("word"):
            lw = w.lower()
            toks.append(("op", _WORD_OPS[lw]) if lw in _WORD_OPS else ("var", w))
        elif const := m.group("const"):
            toks.append(("const", const))
        elif bad := m.group("bad"):
            raise LogicError(f"unexpected character {bad!r} in {expr!r}")
    return toks


class _Parser:
    """Recursive descent; precedence: | then ^ then & then !."""

    def __init__(self, toks: list[tuple[str, str]]):
        self.toks = toks
        self.i = 0
        self.vars: list[str] = []

    def peek(self) -> tuple[str, str] | None:
        return self.toks[self.i] if self.i < len(self.toks) else None

    def eat(self, kind: str, val: str | None = None) -> tuple[str, str]:
        t = self.peek()
        if t is None or t[0] != kind or (val is not None and t[1] != val):
            raise LogicError(f"expected {val or kind}, got {t[1] if t else 'end'}")
        self.i += 1
        return t

    def parse(self):
        node = self.parse_or()
        if self.peek() is not None:
            raise LogicError(f"trailing input at token {self.peek()[1]!r}")
        return node

    def parse_or(self):
        node = self.parse_xor()
        while self.peek() == ("op", "|"):
            self.i += 1
            node = ("or", node, self.parse_xor())
        return node

    def parse_xor(self):
        node = self.parse_and()
        while self.peek() == ("op", "^"):
            self.i += 1
            node = ("xor", node, self.parse_and())
        return node

    def parse_and(self):
        node = self.parse_unary()
        while True:
            t = self.peek()
            if t == ("op", "&"):
                self.i += 1
                node = ("and", node, self.parse_unary())
            elif t and t[0] in {"var", "const", "lparen"}:
                # Implicit AND: "A B" means "A & B".
                node = ("and", node, self.parse_unary())
            else:
                return node

    def parse_unary(self):
        if self.peek() == ("op", "!"):
            self.i += 1
            return ("not", self.parse_unary())
        return self.parse_atom()

    def parse_atom(self):
        t = self.peek()
        if t is None:
            raise LogicError("unexpected end of expression")
        if t[0] == "lparen":
            self.i += 1
            node = self.parse_or()
            self.eat("rparen")
            return node
        if t[0] == "var":
            self.i += 1
            if t[1] not in self.vars:
                self.vars.append(t[1])
            return ("var", t[1])
        if t[0] == "const":
            self.i += 1
            return ("const", t[1] == "1")
        raise LogicError(f"unexpected token {t[1]!r}")


def _eval(node, env: dict[str, bool]) -> bool:
    kind = node[0]
    if kind == "var":
        return env[node[1]]
    if kind == "const":
        return node[1]
    if kind == "not":
        return not _eval(node[1], env)
    a = _eval(node[1], env)
    b = _eval(node[2], env)
    if kind == "and":
        return a and b
    if kind == "or":
        return a or b
    if kind == "xor":
        return a != b
    raise LogicError(f"bad node {kind}")  # pragma: no cover


NAMED_GATES: dict[str, str] = {
    "AND": "A & B",
    "OR": "A | B",
    "NAND": "!(A & B)",
    "NOR": "!(A | B)",
    "XOR": "A ^ B",
    "XNOR": "!(A ^ B)",
    "NOT": "!A",
    "BUFFER": "A",
    "A_NIMPLY_B": "A & !B",
    "A_AND_NOT_B": "A & !B",
    "IMPLY": "!A | B",
    "ALWAYS_ON": "1",
    "ALWAYS_OFF": "0",
    "AND3": "A & B & C",
    "OR3": "A | B | C",
    "AND_NOT": "A & B & !C",
    "MAJORITY3": "(A & B) | (A & C) | (B & C)",
}


@dataclass
class BooleanFunction:
    """One output of a circuit: an ON-set over named inputs."""

    name: str
    inputs: list[str]
    ones: set[int]
    dont_cares: set[int] = frozenset()
    source: str = ""

    @property
    def n(self) -> int:
        return len(self.inputs)

    def evaluate(self, minterm: int) -> bool:
        return minterm in self.ones

    def truth_rows(self) -> list[tuple[str, bool]]:
        return [(format(m, f"0{self.n}b"), m in self.ones) for m in range(1 << self.n)]

    def minimise(self) -> Cover:
        return minimise(self.n, self.ones, self.dont_cares)


def parse_expression(
    expr: str, name: str = "out", inputs: list[str] | None = None
) -> BooleanFunction:
    """Parse a Boolean expression, or a NAMED_GATES key, into a function."""
    key = expr.strip().upper().replace(" ", "_")
    if key in NAMED_GATES:
        expr = NAMED_GATES[key]

    toks = _tokenise(expr)
    if not toks:
        raise LogicError(f"empty expression {expr!r}")
    parser = _Parser(toks)
    tree = parser.parse()

    names = list(inputs) if inputs else parser.vars
    for v in parser.vars:
        if v not in names:
            raise LogicError(f"variable {v!r} not in declared inputs {names}")
    if not names:
        # A constant function still needs one nominal input for a truth table.
        names = ["A"]

    n = len(names)
    ones = {
        m
        for m in range(1 << n)
        if _eval(tree, {nm: bool((m >> (n - 1 - i)) & 1) for i, nm in enumerate(names)})
    }
    return BooleanFunction(name=name, inputs=names, ones=ones, source=expr)


def parse_table(
    spec: str, name: str = "out", inputs: list[str] | None = None
) -> BooleanFunction:
    """Parse a truth table: rows "00:0,01:1", a column "0110", or minterms "1,2"."""
    spec = spec.strip()

    if ":" in spec:
        rows = [r.strip() for r in spec.replace(";", ",").split(",") if r.strip()]
        widths = {len(r.split(":")[0].strip()) for r in rows}
        if len(widths) != 1:
            raise LogicError(f"inconsistent input widths in table {spec!r}")
        n = widths.pop()
        ones, dcs = set(), set()
        for r in rows:
            lhs, _, rhs = r.partition(":")
            lhs = lhs.strip()
            if set(lhs) - set("01"):
                raise LogicError(
                    f"row label {lhs!r} is not binary. Rows are keyed by input bits, "
                    f"so a 2-input table reads '00:0,01:1,10:1,11:0'. To give minterm "
                    f"numbers instead, pass just the list: '1,2'."
                )
            m = int(lhs, 2)
            v = rhs.strip().lower()
            if v in {"1", "on", "true"}:
                ones.add(m)
            elif v in {"x", "-", "dc"}:
                dcs.add(m)
            elif v not in {"0", "off", "false"}:
                raise LogicError(f"bad output value {rhs!r} in row {r!r}")
        names = inputs or _default_names(n)
        return BooleanFunction(name, list(names), ones, frozenset(dcs), spec)

    if set(spec) <= set("01x-"):
        total = len(spec)
        n = max(1, total.bit_length() - 1)
        if (1 << n) != total:
            raise LogicError(
                f"output column has {total} rows, which is not a power of two"
            )
        ones = {m for m, ch in enumerate(spec) if ch == "1"}
        dcs = {m for m, ch in enumerate(spec) if ch in "x-"}
        names = inputs or _default_names(n)
        return BooleanFunction(name, list(names), ones, frozenset(dcs), spec)

    try:
        ms = {int(p) for p in re.split(r"[,\s]+", spec) if p}
    except ValueError as exc:
        raise LogicError(f"cannot parse truth table {spec!r}") from exc
    if not ms:
        raise LogicError(f"empty minterm list {spec!r}")
    n = len(inputs) if inputs else max(1, max(ms).bit_length())
    names = inputs or _default_names(n)
    return BooleanFunction(name, list(names), ms, frozenset(), spec)


def _default_names(n: int) -> list[str]:
    return [chr(ord("A") + i) for i in range(n)]


def all_functions(n: int) -> list[BooleanFunction]:
    """Every Boolean function of n inputs; 16 of them for n = 2."""
    names = _default_names(n)
    return [
        BooleanFunction(
            f"f{mask}",
            list(names),
            {m for m in range(1 << n) if (mask >> m) & 1},
            frozenset(),
            f"mask={mask}",
        )
        for mask in range(1 << (1 << n))
    ]
