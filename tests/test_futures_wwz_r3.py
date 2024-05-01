import subprocess
from os.path import exists

def test_ewkcoffea_r3():
    args = [
        "time",
        "python",
        "analysis/wwz/run_wwz4l.py",
        "-x",
        "futures",
        "input_samples/sample_jsons/test_samples/2022_WWZJetsTo4L2Nu_forCI.json",
        "-o",
        "output_check_yields_2022",
        "-p",
        "analysis/wwz/histos/"
    ]

    # Run ewkcoffea
    subprocess.run(args)

    assert (exists('analysis/wwz/histos/output_check_yields_2022.pkl.gz'))
