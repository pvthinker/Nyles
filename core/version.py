

import nyles
import subprocess


def get_nyles_hash_number():
    """
    Retrieve the hash number of the current nyles version
    """
    gitfolder = nyles.__file__.split('/')[:-2]+['.git']
    with open('/'.join(gitfolder+['HEAD'])) as fid:
        head = fid.readline().split(':')[-1].strip()

    # with open('/'.join(gitfolder+[head])) as fid:
    #     githash = fid.readline().strip()
    return head#githash

def get_git_commit():
    commit = subprocess.check_output("git log -n 1", shell=True)
    status = subprocess.check_output("git status -uno --porcelain", shell=True)
    return commit+status
