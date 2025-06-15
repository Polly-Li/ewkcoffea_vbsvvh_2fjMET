import subprocess
from os.path import exists

### Test for a R2 sample ###

def test_make_yields_after_processor_vvh_r2():
    assert (exists('analysis/vbs_vvh/histos/output_check_yields.pkl.gz')) # Make sure the input pkl file exists

    args = [
        "python",
        "analysis/vbs_vvh/check_vvh_hists.py",
        "analysis/vbs_vvh/histos/output_check_yields.pkl.gz",
        "-o",
        "analysis/vbs_vvh/output_check_yields"
    ]

    # Produce json
    subprocess.run(args)
    assert (exists('analysis/vbs_vvh/output_check_yields.json'))

def test_compare_yields_after_processor_vvh_r2():
    args = [
        "python",
        "ewkcoffea/scripts/comp_json_yields.py",
        "analysis/vbs_vvh/output_check_yields.json",
        "analysis/vbs_vvh/ref_for_ci/counts_vvh_ref.json",
        "-t1",
        "New yields",
        "-t2",
        "Ref yields"
    ]

    # Run comparison
    out = subprocess.run(args, stdout=True)
    assert (out.returncode == 0) # Returns 0 if all pass
