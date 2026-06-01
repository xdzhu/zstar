from __future__ import annotations

from types import MappingProxyType
from typing import Dict, List, Mapping, Tuple

ZERO_TOLERANCE = 1.0e-8
THZ_TO_INV_CM = 33.35641
THZ_TO_MEV = 4.13567
THZ_TO_INV_UM = 1.0e-4 * THZ_TO_INV_CM
DEFAULT_FREQUENCY_UNITS = "?"
DEFAULT_ACTIVITY_UNITS = "AU"
DEFAULT_INTENSITY_UNITS = "AU"

_FREQUENCY_UNIT_LABELS = {
    "thz": "THz",
    "inv_cm": "cm$^{-1}$",
    "mev": "meV",
    "um": r"$\mu$m",
}


def get_frequency_unit_label(frequency_units):
    """Return a TeX-friendly label for the requested frequency unit."""
    if frequency_units is None:
        return DEFAULT_FREQUENCY_UNITS
    return _FREQUENCY_UNIT_LABELS.get(str(frequency_units).lower(), frequency_units)


class _IrrepRecord:
    __slots__ = ("point_group", "ir", "raman", "all")

    def __init__(self, point_group: str, ir: Tuple[str, ...], raman: Tuple[str, ...], all_irreps: Tuple[str, ...]):
        self.point_group = point_group
        self.ir = ir
        self.raman = raman
        self.all = all_irreps

    def as_legacy_dict(self) -> Dict[str, List[str]]:
        return {
            "ir": list(self.ir),
            "raman": list(self.raman),
            "all": list(self.all),
        }


# Format per line:
# point_group ; ir-active irreps ; raman-active irreps ; all irreps
_RAW_ACTIVITY_TABLE = """
1;A;A;A
-1;Au;Ag;Ag,Au
2;A,B;A,B;A,B
m;A',A'';A',A'';A',A''
2/m;Au,Bu;Ag,Bg;Ag,Au,Bg,Bu
222;B1,B2,B3;A,B1,B2,B3;A,B1,B2,B3
mm2;A1,B1,B2;A1,A2,B1,B2;A1,A2,B1,B2
mmm;B1u,B2u,B3u;Ag,B1g,B2g,B3g;Ag,Au,B1g,B1u,B2g,B2u,B3g,B3u
4;A,E;A,B,E;A,B,E
-4;B,E;A,B,E;A,B,E
4/m;Au,Eu;Ag,Bg,Eg;Ag,Au,Bg,Bu,Eg,Eu
422;A2,E;A1,B1,B2,E;A1,A2,B1,B2,E
4mm;A1,E;A1,B1,B2,E;A1,A2,B1,B2,E
-42m;B2,E;A1,B1,B2,E;A1,A2,B1,B2,E
4/mmm;A2u,Eu;A1g,B1g,B2g,Eg;A1g,A1u,A2g,A2u,B1g,B1u,B2g,B2u,Eg,Eu
3;A,E;A,E;A,E
-3;Au,Eu;Ag,Eg;Ag,Au,Eg,Eu
32;A2,E;A1,E;A1,A2,E
3m;A1,E;A1,E;A1,A2,E
-3m;A2u,Eu;A1g,Eg;A1g,A1u,A2g,A2u,Eg,Eu
6;A,E1;A,E1,E2;A,B,E1,E2
-6;A'',E';A',E',E'';A',A'',E',E''
6/m;Au,E1u;Ag,E1g,E2g;Ag,Au,Bg,Bu,E1g,E1u,E2g,E2u
622;A2,E1;A1,E1,E2;A1,A2,B1,B2,E1,E2
6mm;A1,E1;A1,E1,E2;A1,A2,B1,B2,E1,E2
-6m2;A2'',E';A1',E',E'';A1',A1'',A2',A2'',E',E''
6/mmm;A2u,E1u;A1g,E1g,E2g;A1g,A1u,A2g,A2u,B1g,B1u,B2g,B2u,E1g,E1u,E2g,E2u
23;T;A,E,T;A,E,T
m-3;Tu;Ag,Eg,Tg;Ag,Au,Eg,Eu,Tg,Tu
432;T1;A1,E,T2;A1,A2,E,T1,T2
-43m;T2;A1,E,T2;A1,A2,E,T1,T2
m-3m;T1u;A1g,Eg,T2g;A1g,A1u,A2g,A2u,Eg,Eu,T1g,T2g,T1u,T2u
""".strip()


def _split_field(text: str) -> Tuple[str, ...]:
    text = text.strip()
    if not text:
        return tuple()
    return tuple(item.strip() for item in text.split(",") if item.strip())


def _load_registry(table_text: str) -> Dict[str, _IrrepRecord]:
    registry: Dict[str, _IrrepRecord] = {}
    for lineno, raw in enumerate(table_text.splitlines(), start=1):
        line = raw.strip()
        if not line or line.startswith("#"):
            continue
        columns = [part.strip() for part in line.split(";")]
        if len(columns) != 4:
            raise ValueError(f"Malformed activity table at line {lineno}: {raw!r}")
        point_group, ir_text, raman_text, all_text = columns
        point_group = point_group.lower()
        registry[point_group] = _IrrepRecord(
            point_group=point_group,
            ir=_split_field(ir_text),
            raman=_split_field(raman_text),
            all_irreps=_split_field(all_text),
        )
    return registry


_REGISTRY: Mapping[str, _IrrepRecord] = MappingProxyType(_load_registry(_RAW_ACTIVITY_TABLE))
_SUPPORTED_SPECTRA = frozenset(("ir", "raman"))

# Legacy-compatible public variable expected by read_irrep.py fallback branch.
_IRREP_ACTIVITIES: Mapping[str, Dict[str, List[str]]] = MappingProxyType(
    {point_group: record.as_legacy_dict() for point_group, record in _REGISTRY.items()}
)


def _normalize_point_group(point_group) -> str:
    return str(point_group).strip().lower()


def _get_record(point_group) -> _IrrepRecord | None:
    return _REGISTRY.get(_normalize_point_group(point_group))


def _get_inactive_irreps(record: _IrrepRecord, spectrum_key: str) -> List[str]:
    active = set(getattr(record, spectrum_key))
    return [symbol for symbol in record.all if symbol not in active]


def get_irrep_activities(point_group, spectrum_type):
    """
    Return the active and inactive irreps, grouped into 'ir' or 'raman', 
    for the specified point group and spectrum type.

    This preserves the behavior relied on by read_irrep.py.
    """
    normalized_pg = _normalize_point_group(point_group)
    spectrum_key = str(spectrum_type).strip().lower()
    record = _get_record(normalized_pg)

    if record is None:
        print("No activity data for point_group '{0}'.".format(normalized_pg))
        return None

    if spectrum_key not in _SUPPORTED_SPECTRA:
        raise Exception("Error: Unknown spectrum_type '{0}'.".format(spectrum_type))

    active_irreps = list(getattr(record, spectrum_key))
    inactive_irreps = _get_inactive_irreps(record, spectrum_key)
    return active_irreps, inactive_irreps


# Optional helpers. They do not affect existing callers.
def list_supported_point_groups() -> List[str]:
    return list(_REGISTRY.keys())


def has_point_group(point_group) -> bool:
    return _get_record(point_group) is not None


def describe_point_group(point_group) -> Dict[str, List[str]]:
    record = _get_record(point_group)
    if record is None:
        raise KeyError(f"Unknown point group: {point_group}")
    return record.as_legacy_dict()
