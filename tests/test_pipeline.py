from pathlib import Path

from navgen.pipeline import run_generation_pipeline


def test_pipeline_generates_candidates(tmp_path):
    config_path = Path(__file__).resolve().parents[1] / "configs" / "demo.json"
    output_dir = tmp_path / "demo_run"
    result = run_generation_pipeline(config_path, output_dir)

    assert output_dir.joinpath("candidates.csv").exists()
    assert output_dir.joinpath("molecules.sqlite").exists()
    assert len(result["accepted"]) > 0
