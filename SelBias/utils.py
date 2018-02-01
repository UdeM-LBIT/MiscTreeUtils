from hashlib import md5
from signal import signal, SIGINT, SIG_IGN
import os
import numpy as np

def init_worker():
    signal(SIGINT, SIG_IGN)


def load_model(model_obj, tree, path, **kwargs):
    setattr(model_obj, 'run', kwargs.get('run', ''))
    try:
        tree.link_to_evol_model(path, model_obj)
    except KeyError:
        raise(Exception('ERROR: model %s failed, problem with outfile:\n%s' % (
            model_obj.name, path)))


def compute_M(Dn, Ds, inlist, outlist):
    total = []
    for x in inlist:
        for y in outlist:
            if x!=y:
                total.append(Dn[x,y]*1.0/Ds[x,y])
    return np.mean(total)


def name_model(tree, base_name):
    """
    transform the name string into summary of its name and a digestion of the
    full name
    """
    return base_name[:12] + '~' + md5(tree.get_topology_id(attr="name") +
                                      base_name).hexdigest()

def clean_tree(tree):
    """
    remove marks from tree
    """
    for n in tree.get_descendants() + [tree]:
        n.mark = ''


def which(program):
    # keeping it as such
    # too lazy to search for my own "check_binaries"
    def is_exe(fpath):
        return os.path.isfile(fpath) and os.access(fpath, os.X_OK)

    fpath, fname = os.path.split(program)
    if fpath:
        if is_exe(program):
            return program
    else:
        for path in os.environ["PATH"].split(os.pathsep):
            path = path.strip('"')
            exe_file = os.path.join(path, program)
            if is_exe(exe_file):
                return exe_file
    return ""


def find_binary(binary):
    bin_path = which(binary)
    if not os.path.exists(bin_path):
        bin_path = os.path.join(os.path.split(which("ete3"))[0], "ete3_apps", "bin", binary)
    
    if not os.path.exists(bin_path):
        bin_path = os.path.expanduser("~/.etetoolkit/ext_apps-latest/bin/"+binary)

    if not os.path.exists(bin_path):
        raise ValueError("%s binary not found!"%binary)
    return bin_path


class MatrixRep(object):
    def __init__(self, leaflist, defval=0):
            
        self.inlist = leaflist
    
        # keeping tree as key (in case the name was not set for internal nodes)
        self.map = dict((l, i) for i, l in enumerate(self.inlist))
        self.matrix = np.empty((len(self.inlist), len(self.inlist)))
        self.matrix.fill(defval)
        self.shape = self.matrix.shape
    
    def __len__(self):
        """Return the len of the longuest axis"""
        return max(self.shape)
    
    def __contains__(self, item):
        return item in self.inlist
    
    def __iter__(self):
        for (g,i_g) in self.inlist.items():
            for (s,i_s) in self.inlist.items():
                yield (g,s, self.matrix[i_g, i_s])
    
    def _reformat_slice(self, pos, mapping):
        return pos if isinstance(pos, int) else mapping.get(pos, None)

    def _get_new_index(self, index, mapping):
        start =  index.start
        stop =  index.stop
        step = index.step
        start = self._reformat_slice(start, mapping)
        stop = self._reformat_slice(stop, mapping)
        step = step if isinstance(step, int) else None
        return slice(start, stop, step)
        
    def __getitem__(self, index):
        """Indexing with int, string or slice"""
        # the following will return a whole row
        if isinstance(index, basestring):
            return self.matrix[self.map[index]]
        elif isinstance(index, int):
            return self.matrix[index]
        # in th folowing, we are returning a slice
        elif isinstance(index, slice):
            index = self._get_new_index(index, self.map)
            return self.matrix[index] 
        # we are accepting two slices here, no more
        elif len(index) != 2:
            raise TypeError("Invalid index type.")
        # Handle double indexing
        row_index, col_index = index
        if isinstance(row_index, basestring):
            row_index = self.map.get(row_index, None)
        if isinstance(col_index, basestring):
            col_index = self.map.get(col_index, None)        
        elif isinstance(row_index, slice):
            row_index =  self._get_new_index(row_index, self.map)
        elif isinstance(col_index, slice):
            col_index =  self._get_new_index(col_index, self.map)
        # let numpy manage the exceptions
        return self.matrix[row_index, col_index]
    
    def __setitem__(self, index, val):
        """Indexing with int, string or slice"""
        # the following will return a whole row
        if isinstance(index, basestring):
            self.matrix[self.map[index]]  =  val
        elif isinstance(index, int):
            self.matrix[index] = val
        # in th folowing, we are returning a slice
        elif isinstance(index, slice):
            index = self._get_new_index(index, self.map)
            self.matrix[index] = val
        # we are accepting two slices here, no more
        elif len(index) != 2:
            raise TypeError("Invalid index type.")
        # Handle double indexing
        row_index, col_index = index
        if isinstance(row_index, basestring):
            row_index = self.map.get(row_index, None)
        if isinstance(col_index, basestring):
            col_index = self.map.get(col_index, None)        
        elif isinstance(row_index, slice):
            row_index =  self._get_new_index(row_index, self.map)
        elif isinstance(col_index, slice):
            col_index =  self._get_new_index(col_index, self.map)
        # let numpy manage the exceptions
        self.matrix[row_index, col_index] = val