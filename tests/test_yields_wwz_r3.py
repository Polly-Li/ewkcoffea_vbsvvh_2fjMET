import subprocess
from os.path import exists

### Test for a R3 sample ###

def test_make_yields_after_processor_wwz_r3():
    assert (exists('analysis/wwz/histos/output_check_yields_2022.pkl.gz')) # Make sure the input pkl file exists

    args = [
        "python",
        "analysis/wwz/get_wwz_counts.py",
        "-f",
        "analysis/wwz/histos/output_check_yields_2022.pkl.gz",
        "-n",
        "analysis/wwz/output_check_yields_2022",
        "-s",
        "2022_WWZJetsTo4L2Nu",
    ]

    # Produce json
    subprocess.run(args)
    assert (exists('analysis/wwz/output_check_yields.json'))

def test_compare_yields_after_processor_wwz_r3():
    args = [
        "python",
        "ewkcoffea/scripts/comp_json_yields.py",
        "analysis/wwz/output_check_yields_2022.json",
        "analysis/wwz/ref_for_ci/counts_wwz_ref_2022.json",
        "-t1",
        "New yields",
        "-t2",
        "Ref yields"
    ]

    # Run comparison
    out = subprocess.run(args, stdout=True)
    assert (out.returncode == 0) # Returns 0 if all pass
