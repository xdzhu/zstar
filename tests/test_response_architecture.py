import json
from pathlib import Path
import tempfile
import unittest

import numpy as np

from zstar.backends import BackendRegistry, BackendSpec, MetadataBackend, builtin_registry
from zstar.dimensions import DimensionSpec, dimension_spec
from zstar.response_schema import (
    ResponseQuantity,
    ResponseRecord,
    response_record_from_abacus_files,
    response_record_from_bec_result,
    validate_response_document,
)


class ResponseArchitectureTests(unittest.TestCase):
    def test_dimension_defaults_include_physical_one_dimensional_case(self):
        molecule = dimension_spec(0)
        wire = dimension_spec(1)
        slab = dimension_spec(2)
        bulk = dimension_spec(3)
        self.assertEqual(wire.periodic_axes, ("z",))
        self.assertEqual(wire.nonperiodic_axes, ("x", "y"))
        self.assertEqual(wire.intrinsic_response_unit, "angstrom^2")
        self.assertEqual(molecule.intrinsic_response_unit, "angstrom^3")
        self.assertEqual(slab.intrinsic_response_unit, "angstrom")
        self.assertEqual(bulk.intrinsic_response_unit, "1")
        self.assertEqual(dimension_spec(1, "x").periodic_axes, ("x",))

    def test_dimension_axes_must_match_dimension(self):
        with self.assertRaisesRegex(ValueError, "exactly 1"):
            DimensionSpec(1, ("x", "y"))
        with self.assertRaisesRegex(ValueError, "unique"):
            DimensionSpec(2, ("x", "x"))

    def test_response_round_trip_and_validation(self):
        record = ResponseRecord(
            backend="test",
            dimensionality=dimension_spec(1),
            quantities=(
                ResponseQuantity(
                    name="line_polarizability",
                    values=np.eye(3),
                    unit="angstrom^2",
                    normalization="periodic_length",
                    axes=("field", "polarization"),
                ),
            ),
            provenance={"source": "unit-test"},
        )
        with tempfile.TemporaryDirectory() as tmp:
            path = record.write(Path(tmp) / "response.json")
            loaded = ResponseRecord.read(path)
            self.assertEqual(loaded.dimensionality.value, 1)
            self.assertEqual(loaded.quantity("line_polarizability").shape, (3, 3))
            report = validate_response_document(path)
            self.assertTrue(report["valid"])
            self.assertEqual(report["schema_version"], "1.0")

    def test_response_rejects_non_finite_values(self):
        with self.assertRaisesRegex(ValueError, "non-finite"):
            ResponseQuantity(
                name="bad",
                values=[1.0, np.nan],
                unit="1",
                normalization="none",
            )

    def test_existing_bec_result_can_be_normalized(self):
        source = {
            "schema_version": 1,
            "backend": "vasp",
            "method": "dfpt",
            "tensor_convention": "rows=displacement; columns=polarization",
            "epsilon_infinity": np.eye(3).tolist(),
            "atoms": [
                {"index": 1, "label": "Si", "tensor": (np.eye(3) * 2).tolist()},
                {"index": 2, "label": "C", "tensor": (np.eye(3) * -2).tolist()},
            ],
            "acoustic_sum_tensor": np.zeros((3, 3)).tolist(),
        }
        record = response_record_from_bec_result(
            source,
            dimensionality=2,
            provenance={"file": "vasp_bec.json"},
        )
        self.assertEqual(record.dimensionality.periodic_axes, ("x", "y"))
        self.assertEqual(record.quantity("born_effective_charge").shape, (2, 3, 3))
        self.assertEqual(
            record.quantity("supercell_electronic_dielectric").metadata[
                "intrinsic_low_dimensional_response_required"
            ],
            True,
        )

    def test_registry_reports_only_implemented_capabilities(self):
        registry = builtin_registry()
        self.assertTrue(registry.get("pyatb").spec.supports("born_effective_charge", 2))
        self.assertFalse(registry.get("vasp").spec.supports("born_effective_charge", 2))
        self.assertTrue(registry.get("phonopy").spec.supports("gamma_modes", 1))
        custom = BackendRegistry()
        custom.register(
            MetadataBackend(
                BackendSpec(
                    name="demo",
                    display_name="Demo",
                    capabilities={"forces": frozenset({0, 1, 2, 3})},
                )
            )
        )
        self.assertEqual(custom.names(), ("demo",))

    def test_abacus_legacy_products_can_be_normalized(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            zborn = root / "Z-BORN-all.out"
            zborn.write_text(
                "No. Atom xx xy xz yx yy yz zx zy zz\n"
                "1 Si 2 0 0 0 2 0 0 0 2\n",
                encoding="utf-8",
            )
            born = root / "BORN"
            born.write_text(
                "5 0 0 0 5 0 0 0 5\n2 0 0 0 2 0 0 0 2\n",
                encoding="utf-8",
            )
            record = response_record_from_abacus_files(
                zborn, born_path=born, dimensionality=2
            )
            self.assertEqual(record.backend, "abacus")
            self.assertEqual(record.dimensionality.value, 2)
            self.assertEqual(record.quantity("born_effective_charge").shape, (1, 3, 3))
            self.assertEqual(
                record.quantity("supercell_electronic_dielectric").shape, (3, 3)
            )


if __name__ == "__main__":
    unittest.main()
