import os
import ewkcoffea

pjoin = os.path.join

# This function takes as input any path (inside of ewkcoffea/ ewkcoffea), and returns the absolute path
def ewkcoffea_path(path_in_repo):
    return pjoin(ewkcoffea.__path__[0], path_in_repo)
