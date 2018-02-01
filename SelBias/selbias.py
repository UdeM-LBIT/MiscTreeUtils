# Compute dn/ds and mx/mo (see Nadia's paper) on a tree
# ~30% of this code is adapted from the ete3 evol module

from ete3 import EvolTree
from ete3.evol.control import AVAIL
from ete3.evol import Model
from multiprocessing import Pool, Queue
from subprocess import Popen, PIPE
from utils import *
import logging
import os
import collections
import argparse
import uuid


class ValidateModels(argparse.Action):

    def __call__(self, parser, args, values, option_string=None):
        valid_test = [('bsA1', 'bsA'), ('M0', 'b_free'), ('b_neut', 'b_free'),
                      ('M3', 'bsD'), ('M1', 'bsA')]
        # (bsa1 vs bsa) is the default positive selection
        h0, h1 = values.split('|')
        if (h0, h1) not in valid_test:
            raise ValueError('Invalid Hypothesis testing %s vs %s' % (h0, h1))
        setattr(args, self.dest, [h0, h1])


def check_done(tree, modmodel, results):
    dir_name = name_model(tree, modmodel)
    if os.path.exists(os.path.join(tree.workdir, dir_name, 'out')):
        fhandler = open(os.path.join(tree.workdir, dir_name, 'out'))
        fhandler.seek(-50, 2)
        if 'Time used' in fhandler.read():
            model_obj = Model(modmodel, tree)

            results.append((os.path.join(tree.workdir, dir_name, "out"),
                            model_obj, ''))
            return True
    return False


def run_models(tree, models, maxproc, binary, freemod=False, **kwargs):
    # TO BE IMPROVED: multiprocessing should be called in a simpler way
    print("\nRunning CodeML (%s CPUs)" % maxproc)
    pool = Pool(maxproc or None, init_worker)
    binary = os.path.expanduser(binary)

    results = []
    for model in models:
        logging.debug('  - processing model %s' % model)
        if AVAIL[model.split('.')[0]]['allow_mark']:
            for node in tree.iter_descendants():
                # only on internal nodes ?
                branch = [n.node_id for n in node.traverse()]
                clean_tree(tree)
                tree.mark_tree(branch, marks=['#1']*len(branch))
                modmodel = model + '.' + str(node.node_id)
                if check_done(tree, modmodel, results):
                    logging.debug('Model %s (%s) already executed... SKIPPING' %
                          (model,  name_model(tree, modmodel)))
                    continue
                results.append(pool.apply_async(
                    local_run_model, args=(tree.copy(), modmodel, binary), kwds=kwargs))
        else:
            if check_done(tree, model, results):
                logging.debug('Model %s (%s) already executed... SKIPPING' % (
                    model, name_model(tree, model)))
                continue
            results.append(pool.apply_async(
                local_run_model, args=(tree.copy(), model, binary), kwds=kwargs))

    pool.close()
    pool.join()

    clean_tree(tree)
    # join back results to tree
    for result in results:
        try:
            path, model_obj, info = result.get()
        except AttributeError:
            path, model_obj, info = result
        load_model(model_obj, tree, path, run=info, **kwargs)

    for node in tree.iter_descendants():
        if freemod:
            node.add_feature("dnds", tree.get_evol_model(
                models[0]).branches[node.node_id].get('w', None))

        elif len(models) == 3:
            M0, bneut, bfree = ['M0']+[x+"."+str(node.node_id) for x in models[1:]]
            logging.debug('\n==> Node: %s'%node.node_id)

            logging.debug(node.get_ascii(show_internal=True))
            logging.debug('\n---------------M0 Model---------------------\n')
            logging.debug(tree.get_evol_model(M0))
            # and also more info on negative selection test
            logging.debug('\n---------------b_neut Model---------------------\n')
            logging.debug(tree.get_evol_model(bneut))
            logging.debug('\n---------------b_free Model---------------------\n')
            logging.debug(tree.get_evol_model(bfree))
            logging.debug('\n')

            has_diff_ratio = tree.get_most_likely(bfree, M0)
            best_model = M0
            message = "One ratio for all branches"
            if has_diff_ratio < 0.05:
                # we have selection probably
                best_model = bfree
                under_selection = tree.get_most_likely(bfreee, bneut)
                message = "different ratio on foreground branch (relaxation)"
                if under_selection >= 0.05:
                    # this is more like relaxation
                    best_model = bneut
                else:
                    message = "Possible selection on foreground branch"
            
            logging.debug('xxxxxxx %s xxxxxx\n'%message)

            omega = {
                M0: tree.get_evol_model(M0).branches[node.node_id]['w'],
                bneut: tree.get_evol_model(bneut).branches[node.node_id]['w'],
                bfree: tree.get_evol_model(bfree).branches[node.node_id]['w'],
            }
            node.add_feature("best", best_model.split('.')[0])
            node.add_feature("omega", omega)
            node.add_feature("dnds", omega[best_model])

        else:
            # this is hypothesis testing
            null, altn = models
            null += '.' + str(node.node_id)
            altn += '.' + str(node.node_id)
            res = tree.get_most_likely(altn, null)
            best_model = (null if res >= 0.05 else altn)
            
            logging.debug('\n==> Node: %s'%node.node_id)
            logging.debug(node.get_ascii(show_internal=True))
            logging.debug('\n---------------Null Model---------------------\n')
            logging.debug(tree.get_evol_model(null))
            # and also more info on negative selection test
            logging.debug('\n---------------Alt Model---------------------\n')
            logging.debug(tree.get_evol_model(altn))
            logging.debug('\n')
            message = "failed to reject null model (%s)"%null
            if best_model == altn:
                message = "null model (%s) rejected"%null
            logging.debug('xxxxxxx %s xxxxxx\n'%message)

            null_pos, max_null = max(enumerate(tree.get_evol_model(null).classes['proportions']), key=lambda x: x[1])
            altn_pos, max_altn = max(enumerate(tree.get_evol_model(altn).classes['proportions']), key=lambda x: x[1])

            omega = {
                null: tree.get_evol_model(null).classes['foreground w'][null_pos],
                altn: tree.get_evol_model(altn).classes['foreground w'][altn_pos]
            }
            node.add_feature("best", best_model.split('.')[0])
            node.add_feature("omega", omega)
            node.add_feature("dnds", omega[best_model])

    return tree


def local_run_model(tree, model_name, binary, mxmo=False, model_obj=None, **kwargs):
    '''
    local version of model runner. Needed for multiprocessing pickling...
    '''
    def clean_exit(a, b):
        if proc:
            logging.debug("Killing process %s" % proc)
            proc.terminate()
            proc.kill(-9)
        exit(a, b)
    proc = None
    signal(SIGINT, clean_exit)

    if not model_obj:
        model_obj = Model(model_name, tree, **kwargs)
    # dir_name = model_obj.name
    fullpath = os.path.join(tree.workdir, name_model(tree, model_obj.name))
    os.system("mkdir -p %s" % fullpath)
    # write algn file
    tree._write_algn(fullpath + '/algn')
    # write tree file
    tree.write(outfile=fullpath+'/tree',
               format=(10 if model_obj.properties['allow_mark'] else 9))
    ctrl_string = model_obj.get_ctrl_string(fullpath+'/tmp.ctl')

    hlddir = os.getcwd()

    os.chdir(fullpath)

    proc = Popen("%s tmp.ctl" % binary, stdout=PIPE, stdin=PIPE, shell=True)
    proc.stdin.write('\n')  # in case codeml/slr asks something
    job, err = proc.communicate()
    if err is not None or b'error' in job or b'Error' in job:
        print("ERROR: inside CodeML!!\n" + job)
        return (None, None)
    os.chdir(hlddir)

    if mxmo:
        return fullpath, None, job
    return os.path.join(fullpath, 'out'), model_obj, job


def mat_parse(infile, tlist):
    mat = MatrixRep(tlist)
    with open(infile) as IN:
        first = True
        max_item = 0
        namelist = []
        for line in IN:
            line = line.strip()
            if line:
                if first:
                    max_item = int(line)
                    first = False
                else:
                    line = line.split()
                    name = line[0]
                    namelist.append(name)
                    for pos, val in enumerate(line[1:]):
                        mat[name, namelist[pos]] = val
                        mat[namelist[pos], name] = val

    return mat


def run_mxmo(tree, codeml_bin):
    mod = Model('M0.pairwise', tree, runmode=-2, cleandata=1)
    fullpath, _, info = local_run_model(
        tree, mod.name, codeml_bin, mxmo=True, model_obj=mod)
    leafset = tree.get_leaf_names()
    Dn = mat_parse(os.path.join(fullpath, '2ML.dN'), leafset)
    Ds = mat_parse(os.path.join(fullpath, '2ML.dS'), leafset)
    for node in tree.iter_descendants("postorder"):
        if not node.is_leaf():
            innode = node.get_leaf_names()
            outnode = (set(leafset) - set(innode))
            Mo = compute_M(Dn, Ds, innode, outnode)
            Mx = compute_M(Dn, Ds, innode, innode)
            node.add_feature('mxmo', Mo*1.0/Mx)
    return tree


def run(tree, freemod=False, mxmo=True, codeml=None, hyptest=None, cores=0, **kwargs):
    # in case we only got 1 model :(

    codeml_bin = find_binary(codeml)
    if mxmo:
        tree = run_mxmo(tree, codeml_bin)
    else:
        models = []
        if freemod:
            models = ['fb']
        elif hyptest:
            models = hyptest
            if ['b_neut', 'b_free'] == list(models):
                models = ['M0'] + ['b_neut', 'b_free']

        tree = run_models(tree, models, cores, codeml_bin,
                          freemod=freemod, **kwargs)

    outfile = os.path.join(tree.workdir, 'output')
    with open(outfile, 'w') as OUT:
        OUT.write('## tree\n')
        OUT.write(tree.write(format=0, features=['node_id']))
        OUT.write('\n## info\n')
        OUT.write("\t".join(['id', 'mx/mo', 'dn/ds', 'bestmod'])+"\n")

        # lazyness overload
        for node in tree.traverse():
            if not hasattr(node, 'mxmo'):
                node.add_features(mxmo='-')
            else:
                node.mxmo = "%.3f"%node.mxmo
            if not hasattr(node, 'dnds'):
                node.add_features(dnds='-')
            else:
                node.dnds = "%.3f"%node.dnds
            if not hasattr(node, 'best'):
                node.add_features(best='-')
            OUT.write("\t".join([str(node.node_id), node.mxmo, node.dnds, str(node.best)])+"\n")

    print("Results saved in : %s"%outfile)

def exec_program(cmdline):
    proc = Popen(cmdline, stdout=PIPE, stdin=PIPE, shell=True)
    proc.stdin.write('\n')
    job, err = proc.communicate()
    if err is not None or b'error' in job or b'Error' in job:
        print("ERROR while running %s!!\n%s" % (cmdline, job))
        return False
    return True


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='SelBias')
    parser.add_argument('-t', '--gtree', dest='genetree', required=True,
                        help="Either the filename or the newick of the genetree")
    parser.add_argument('-w', '--wdir', dest="wdir",
                        help="Working directory to use. Default will be a tmp folder created in the current directory")
    parser.add_argument('--gcode', default='1', help="Genetic code to use")

    algo = parser.add_mutually_exclusive_group(required=True)    
    algo.add_argument('--mxmo',  action="store_true",
                        help="Compute MxMo as described in the paper")
    algo.add_argument('--hyptest', metavar='H0|H1', action=ValidateModels,
                        help="The two models to test. "
                        "List of valid pair are:\n "
                        "bsA1|bsA, M0|b_free, b_neut|b_free, M3|bsD, M1|bsA. "
                        "Just use 'b_neut|b_free') or avoid this option.")
    algo.add_argument('--freemodel', action="store_true",
                        help="Free ratio on each branch of the tree")

    seqtype = parser.add_mutually_exclusive_group(required=True)
    seqtype.add_argument('--codaln', '-c', dest='codaln',
                         help="Codon alignment")
    seqtype.add_argument('--seq', '-s', dest='seq', nargs=2, metavar=(
        'protaln', 'nucseq'), help="nucleotide sequence and protein alignment as an alternative to codon alignment")

    parser.add_argument('--codeml', default="codeml",
                        help="Location of the codeml binary. It is easier to add it to your home path")
    parser.add_argument('--debug', action="store_true", help="Print debug information")
    parser.add_argument('--alnformat', '-f', dest="alnformat", default="paml",
                        help="Alignment format to read. Default is paml, which is the format returned by pal2nal")
    parser.add_argument("--pal2nal_args", dest="pal2nal",
                        default="", help="Supplementary arguments for pal2nal")
    parser.add_argument("--cpu", dest="cores", type=int,
                        default=1, help="Maximum number of CPU cores"
                        " available in the execution host. If higher"
                        " than 1, multi-threading will be enabled"
                        " (if 0 all available cores will be used)")

    args = parser.parse_args()
    session_id = uuid.uuid4().hex[:6]
    if args.debug:
        logging.basicConfig(format='%(message)s', level=logging.DEBUG)

    if not (args.mxmo or args.freemodel or args.hyptest):
        raise ValueError(
            "You have to provide one way for computing dn/ds: [ --mxmo | --hyptest | --freemodel")

    if not args.wdir:
        print("Current session id is %s" % session_id)
        args.wdir = os.path.join(os.getcwd(), 'selbias.res', session_id)

    if not os.path.exists(args.wdir):
        os.makedirs(args.wdir)

    codaln = None
    if args.seq:
        codaln = os.path.join(args.wdir, 'cod.aln')
        # build codon alignment from aa alignment and nuc sequence
        pal2nal_cmd = "bin/pal2nal.pl %s %s -codontable %s -output paml %s > %s" % (
            args.seq[0], args.seq[1], args.gcode, args.pal2nal, codaln)
        print pal2nal_cmd
        success = exec_program(pal2nal_cmd)
        if not success:
            raise ValueError("Failed to run pal2nal to convert alignment")
        args.alnformat = 'paml'
    else:
        codaln = args.codaln

    tree = EvolTree(args.genetree)
    tree.link_to_alignment(codaln, alg_format=args.alnformat)
    tree.workdir = args.wdir

    run(tree, freemod=args.freemodel, mxmo=args.mxmo,
        codeml=args.codeml, hyptest=args.hyptest, cores=args.cores)
