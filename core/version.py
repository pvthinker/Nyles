

import nyles


def get_nyles_hash_number():
    """
    Retrieve the hash number of the current nyles version
    """
    gitfolder = nyles.__file__.split('/')[:-2]+['.git']
    with open('/'.join(gitfolder+['HEAD'])) as fid:
        head = fid.readline().split(':')[-1].strip()

    with open('/'.join(gitfolder+[head])) as fid:
        githash = fid.readline().strip()
    return githash
