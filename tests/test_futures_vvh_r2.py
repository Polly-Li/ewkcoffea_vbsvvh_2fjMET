import subprocess
from os.path import exists

def test_ewkcoffea_vvh():
    args = [
        "time",
        "python",
        "analysis/vbs_vvh/run_analysis.py",
        "-x",
        "futures",
        "input_samples/sample_jsons/test_samples/UL17_WWZJetsTo4L2Nu_forCI.json",
        "-o",
        "output_check_yields",
        "-outpath",
        "analysis/vbs_vvh/histos/"
        "-p",
        "1lep1jf_nano",
    ]

    # Run ewkcoffea
    subprocess.run(args)

    assert (exists('analysis/vbs_vvh/histos/output_check_yields.pkl.gz'))
