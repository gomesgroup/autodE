import hashlib
import json

import pytest
from scipy.optimize import OptimizeResult

from autode.atoms import Atom
from autode.neb.ci import CINEB
from autode.species import Molecule
from autode.wrappers.ORCA import ORCA
from autode.wrappers.XTB import XTB
from autode.wrappers.mlip_external import (
    MLIPAcceleratedNEB,
    MLIPNEBResult,
)


@pytest.fixture
def endpoints():
    reactant = Molecule(
        name="reactant",
        atoms=[Atom("H"), Atom("H", x=0.8)],
        mult=1,
    )
    product = Molecule(
        name="product",
        atoms=[Atom("H"), Atom("H", x=1.2)],
        mult=1,
    )
    return reactant, product


def test_mlip_neb_requires_a_concrete_method(endpoints, tmp_path):
    reactant, product = endpoints

    with pytest.raises(TypeError, match="Method"):
        MLIPAcceleratedNEB(
            reactant,
            product,
            method="uma-small",
            method_provenance={"model_identity": "UMA Small v1.2"},
            result_dir=tmp_path,
        )


def test_mlip_neb_requires_nonempty_serializable_provenance(
    endpoints, tmp_path
):
    reactant, product = endpoints
    method = XTB()

    with pytest.raises(ValueError, match="provenance"):
        MLIPAcceleratedNEB(
            reactant,
            product,
            method=method,
            method_provenance={},
            result_dir=tmp_path,
        )

    with pytest.raises(ValueError, match="JSON"):
        MLIPAcceleratedNEB(
            reactant,
            product,
            method=method,
            method_provenance={"callback": lambda: None},
            result_dir=tmp_path,
        )


def test_mlip_neb_requires_at_least_three_images(endpoints, tmp_path):
    reactant, product = endpoints

    with pytest.raises(ValueError, match="three"):
        MLIPAcceleratedNEB(
            reactant,
            product,
            method=XTB(),
            method_provenance={"model_identity": "UMA Small v1.2"},
            n_images=2,
            result_dir=tmp_path,
        )


def test_mlip_neb_stores_explicit_method_and_provenance(endpoints, tmp_path):
    reactant, product = endpoints
    method = XTB()
    provenance = {
        "model_alias": "uma-s-1p2",
        "model_identity": "UMA Small v1.2",
        "progext_sha256": "abc123",
    }

    neb = MLIPAcceleratedNEB(
        reactant,
        product,
        method=method,
        method_provenance=provenance,
        n_images=7,
        result_dir=tmp_path,
    )

    assert neb.method is method
    assert neb.method_provenance == provenance
    assert neb.n_images == 7
    assert neb.result_dir == tmp_path.resolve()


def test_dft_refinement_is_unconditionally_fail_closed():
    neb = object.__new__(MLIPAcceleratedNEB)
    neb.ts_guess = object()

    with pytest.raises(NotImplementedError, match="MLIP-only"):
        neb.refine_with_dft()


def test_result_type_is_exported_from_wrappers_namespace():
    import autode.wrappers as wrappers

    assert wrappers.MLIPNEBResult is MLIPNEBResult


class FakeCINEB:
    def __init__(
        self,
        peak_species,
        optimize_result=None,
        error=None,
        images=None,
    ):
        self.peak_species = peak_species
        self.optimize_result = optimize_result
        self.error = error
        self.calculate_kwargs = None
        self.images = [] if images is None else images

    def calculate(self, **kwargs):
        self.calculate_kwargs = kwargs
        if self.error is not None:
            raise self.error
        return self.optimize_result


def _configured_neb(endpoints, tmp_path, n_images=5):
    reactant, product = endpoints
    return MLIPAcceleratedNEB(
        reactant,
        product,
        method=XTB(),
        method_provenance={
            "model_alias": "uma-s-1p2",
            "model_identity": "UMA Small v1.2",
        },
        n_images=n_images,
        result_dir=tmp_path,
    )


def _replace_cineb_constructor(monkeypatch, expected_endpoints, fake_cineb):
    expected_reactant, expected_product = expected_endpoints

    def from_end_points(cls, initial, final, num):
        assert initial is expected_reactant
        assert final is expected_product
        assert num == 5
        return fake_cineb

    monkeypatch.setattr(CINEB, "from_end_points", classmethod(from_end_points))


def _checkpoint_images():
    images = [
        Molecule(
            name=f"image_{index}",
            atoms=[Atom("H"), Atom("H", x=distance)],
            mult=1,
        )
        for index, distance in enumerate((0.8, 0.9, 1.0, 1.1, 1.2))
    ]
    for image, energy in zip(images, (0.0, 0.2, 0.4, 0.3, 0.1)):
        image.energy = energy
    return images


def test_run_mlip_neb_uses_native_cineb_and_explicit_runner(
    endpoints, tmp_path, monkeypatch
):
    peak = Molecule(
        name="path_peak",
        atoms=[Atom("H"), Atom("H", x=1.0)],
        mult=1,
    )
    peak.energy = 0.4
    fake_cineb = FakeCINEB(
        peak_species=peak,
        optimize_result=OptimizeResult(
            success=True,
            message="converged",
            nit=7,
            nfev=19,
        ),
        images=_checkpoint_images(),
    )
    _replace_cineb_constructor(monkeypatch, endpoints, fake_cineb)
    neb = _configured_neb(endpoints, tmp_path)
    calculation_runner = object()

    result = neb.run_mlip_neb(
        max_steps=123,
        n_cores=4,
        calculation_runner=calculation_runner,
        resume=False,
    )

    assert fake_cineb.calculate_kwargs == {
        "method": neb.method,
        "n_cores": 4,
        "max_n": 123,
        "calculation_runner": calculation_runner,
    }
    assert neb.mlip_path is fake_cineb
    assert result.status == "succeeded"
    assert result.converged is True
    assert result.optimizer_message == "converged"
    assert result.optimizer_n_iterations == 7
    assert result.optimizer_n_evaluations == 19
    assert result.n_images == 5
    assert result.n_cores == 4
    assert result.max_steps == 123
    assert result.method_provenance["model_identity"] == "UMA Small v1.2"
    assert neb.get_ts_guess() is not None
    assert neb.get_ts_guess().coordinates.tolist() == peak.coordinates.tolist()


@pytest.mark.parametrize(
    ("optimize_result", "peak_species", "expected_status"),
    [
        (
            OptimizeResult(
                success=False,
                message="maximum evaluations reached",
                nit=10,
                nfev=20,
            ),
            Molecule(atoms=[Atom("H"), Atom("H", x=1.0)], mult=1),
            "nonconverged",
        ),
        (
            OptimizeResult(
                success=True,
                message="converged",
                nit=4,
                nfev=9,
            ),
            None,
            "no_peak",
        ),
    ],
)
def test_run_mlip_neb_returns_clear_unsuccessful_state(
    endpoints,
    tmp_path,
    monkeypatch,
    optimize_result,
    peak_species,
    expected_status,
):
    fake_cineb = FakeCINEB(
        peak_species=peak_species,
        optimize_result=optimize_result,
        images=_checkpoint_images(),
    )
    _replace_cineb_constructor(monkeypatch, endpoints, fake_cineb)
    neb = _configured_neb(endpoints, tmp_path)

    result = neb.run_mlip_neb(resume=False)

    assert result.status == expected_status
    assert result.converged is bool(optimize_result.success)
    assert result.ts_guess_coordinates is None
    assert neb.get_ts_guess() is None


def test_run_mlip_neb_records_failed_state_then_reraises(
    endpoints, tmp_path, monkeypatch
):
    fake_cineb = FakeCINEB(
        peak_species=None,
        error=RuntimeError("gateway evaluation failed"),
        images=_checkpoint_images(),
    )
    _replace_cineb_constructor(monkeypatch, endpoints, fake_cineb)
    neb = _configured_neb(endpoints, tmp_path)

    with pytest.raises(RuntimeError, match="gateway evaluation failed"):
        neb.run_mlip_neb(resume=False)

    assert neb.last_result.status == "failed"
    assert neb.last_result.error_type == "RuntimeError"
    assert neb.last_result.error_message == "gateway evaluation failed"
    assert neb.get_ts_guess() is None


def test_path_construction_failure_persists_clear_failed_state(
    endpoints, tmp_path, monkeypatch
):
    monkeypatch.setattr(
        CINEB,
        "from_end_points",
        classmethod(
            lambda cls, *args, **kwargs: (_ for _ in ()).throw(
                RuntimeError("IDPP construction failed")
            )
        ),
    )
    neb = _configured_neb(endpoints, tmp_path)

    with pytest.raises(RuntimeError, match="IDPP construction failed"):
        neb.run_mlip_neb(resume=False)

    receipt = json.loads(
        (tmp_path / "mlip-neb-result.json").read_text()
    )
    assert receipt["status"] == "failed"
    assert receipt["error_type"] == "RuntimeError"
    assert receipt["checkpoint_path"] is None
    assert neb.last_result.status == "failed"


@pytest.mark.parametrize(
    ("kwargs", "message"),
    [
        ({"max_steps": 0}, "max_steps"),
        ({"n_cores": 0}, "n_cores"),
    ],
)
def test_run_mlip_neb_rejects_nonpositive_controls(
    endpoints, tmp_path, kwargs, message
):
    neb = _configured_neb(endpoints, tmp_path)

    with pytest.raises(ValueError, match=message):
        neb.run_mlip_neb(resume=False, **kwargs)


def test_run_persists_atomic_provenance_and_checkpoint(
    endpoints, tmp_path, monkeypatch
):
    images = _checkpoint_images()
    fake_cineb = FakeCINEB(
        peak_species=images[1],
        optimize_result=OptimizeResult(
            success=True,
            message="converged",
            nit=5,
            nfev=11,
        ),
        images=images,
    )
    _replace_cineb_constructor(monkeypatch, endpoints, fake_cineb)
    neb = _configured_neb(endpoints, tmp_path)

    result = neb.run_mlip_neb(
        max_steps=41,
        n_cores=4,
        calculation_runner=object(),
        resume=False,
    )

    receipt_path = tmp_path / "mlip-neb-result.json"
    checkpoint_path = tmp_path / "mlip-neb-checkpoint.xyz"
    assert receipt_path.is_file()
    assert checkpoint_path.is_file()
    receipt = json.loads(receipt_path.read_text())
    assert receipt["schema_version"] == 1
    assert receipt["status"] == "succeeded"
    assert receipt["request"]["n_images"] == 5
    assert receipt["request"]["max_steps"] == 41
    assert receipt["request"]["n_cores"] == 4
    assert receipt["request"]["method_provenance"] == (
        neb.method_provenance
    )
    assert receipt["checkpoint_sha256"] == hashlib.sha256(
        checkpoint_path.read_bytes()
    ).hexdigest()
    assert result.checkpoint_path == str(checkpoint_path.resolve())
    assert result.checkpoint_sha256 == receipt["checkpoint_sha256"]
    assert not list(tmp_path.glob(".*.tmp.*"))


def test_successful_matching_result_resumes_without_recalculation(
    endpoints, tmp_path, monkeypatch
):
    images = _checkpoint_images()
    first_path = FakeCINEB(
        peak_species=images[1],
        optimize_result=OptimizeResult(
            success=True,
            message="converged",
            nit=5,
            nfev=11,
        ),
        images=images,
    )
    _replace_cineb_constructor(monkeypatch, endpoints, first_path)
    _configured_neb(endpoints, tmp_path).run_mlip_neb(
        max_steps=41,
        n_cores=4,
        resume=False,
    )

    resumed_path = FakeCINEB(
        peak_species=images[1],
        optimize_result=None,
        images=images,
    )

    def from_file(cls, filename):
        assert filename == str(
            (tmp_path / "mlip-neb-checkpoint.xyz").resolve()
        )
        return resumed_path

    monkeypatch.setattr(CINEB, "from_file", classmethod(from_file))
    monkeypatch.setattr(
        CINEB,
        "from_end_points",
        classmethod(
            lambda cls, *args, **kwargs: pytest.fail(
                "matching resume reconstructed endpoints"
            )
        ),
    )
    neb = _configured_neb(endpoints, tmp_path)

    result = neb.run_mlip_neb(
        max_steps=41,
        n_cores=4,
        resume=True,
    )

    assert result.status == "succeeded"
    assert result.resumed is True
    assert resumed_path.calculate_kwargs is None
    assert neb.mlip_path is resumed_path
    assert neb.get_ts_guess() is not None


def test_nonconverged_checkpoint_resumes_native_cineb(
    endpoints, tmp_path, monkeypatch
):
    images = _checkpoint_images()
    first_path = FakeCINEB(
        peak_species=images[1],
        optimize_result=OptimizeResult(
            success=False,
            message="maximum evaluations reached",
            nit=8,
            nfev=17,
        ),
        images=images,
    )
    _replace_cineb_constructor(monkeypatch, endpoints, first_path)
    first = _configured_neb(endpoints, tmp_path)
    assert first.run_mlip_neb(resume=False).status == "nonconverged"

    resumed_path = FakeCINEB(
        peak_species=images[1],
        optimize_result=OptimizeResult(
            success=True,
            message="converged after resume",
            nit=3,
            nfev=7,
        ),
        images=images,
    )
    monkeypatch.setattr(
        CINEB,
        "from_file",
        classmethod(lambda cls, filename: resumed_path),
    )
    monkeypatch.setattr(
        CINEB,
        "from_end_points",
        classmethod(
            lambda cls, *args, **kwargs: pytest.fail(
                "nonconverged resume reconstructed endpoints"
            )
        ),
    )
    neb = _configured_neb(endpoints, tmp_path)
    runner = object()

    result = neb.run_mlip_neb(
        calculation_runner=runner,
        resume=True,
    )

    assert result.status == "succeeded"
    assert result.resumed is True
    assert resumed_path.calculate_kwargs["calculation_runner"] is runner


@pytest.mark.parametrize(
    "mutation",
    ["provenance", "endpoint", "method"],
)
def test_resume_rejects_request_identity_mismatch(
    endpoints, tmp_path, monkeypatch, mutation
):
    images = _checkpoint_images()
    first_path = FakeCINEB(
        peak_species=images[1],
        optimize_result=OptimizeResult(
            success=True,
            message="converged",
            nit=5,
            nfev=11,
        ),
        images=images,
    )
    _replace_cineb_constructor(monkeypatch, endpoints, first_path)
    _configured_neb(endpoints, tmp_path).run_mlip_neb(resume=False)

    reactant, product = endpoints
    method = XTB()
    provenance = {
        "model_alias": "uma-s-1p2",
        "model_identity": "UMA Small v1.2",
    }
    if mutation == "provenance":
        provenance["checkpoint_sha256"] = "changed"
    elif mutation == "endpoint":
        reactant = reactant.copy()
        reactant.atoms[0].translate([0.1, 0.0, 0.0])
    else:
        method = ORCA()
    neb = MLIPAcceleratedNEB(
        reactant,
        product,
        method=method,
        method_provenance=provenance,
        n_images=5,
        result_dir=tmp_path,
    )

    with pytest.raises(ValueError, match="resume request"):
        neb.run_mlip_neb(resume=True)


def test_resume_rejects_checkpoint_digest_mismatch(
    endpoints, tmp_path, monkeypatch
):
    images = _checkpoint_images()
    first_path = FakeCINEB(
        peak_species=images[1],
        optimize_result=OptimizeResult(
            success=True,
            message="converged",
            nit=5,
            nfev=11,
        ),
        images=images,
    )
    _replace_cineb_constructor(monkeypatch, endpoints, first_path)
    _configured_neb(endpoints, tmp_path).run_mlip_neb(resume=False)
    checkpoint = tmp_path / "mlip-neb-checkpoint.xyz"
    checkpoint.write_text(checkpoint.read_text() + "\ncorrupted\n")

    with pytest.raises(ValueError, match="checkpoint digest"):
        _configured_neb(endpoints, tmp_path).run_mlip_neb(resume=True)


def test_checkpoint_round_trips_through_native_cineb(
    endpoints, tmp_path
):
    neb = _configured_neb(endpoints, tmp_path)
    images = _checkpoint_images()
    images[2].atoms[1].coord[0] = 1.012345678901234
    neb.mlip_path = CINEB.from_list(images)

    checkpoint_path, checkpoint_sha256 = (
        neb._atomic_write_checkpoint()
    )
    loaded = CINEB.from_file(checkpoint_path)

    assert len(loaded.images) == 5
    assert [float(image.energy) for image in loaded.images] == pytest.approx(
        [0.0, 0.2, 0.4, 0.3, 0.1]
    )
    assert loaded.images[2].atoms[1].coord[0] == pytest.approx(
        1.012345678901234,
        abs=1e-12,
    )
    assert checkpoint_sha256 == hashlib.sha256(
        (tmp_path / "mlip-neb-checkpoint.xyz").read_bytes()
    ).hexdigest()
