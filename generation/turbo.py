#! /bin/sh
""":" .

exec python "$0" ${1+"$@"}
"""

import atexit
import datetime
import os
import pickle
import shutil
import subprocess
import sys
import tempfile

import pathforcer

__doc__ = """Generate rules that reduce the complexity."""

# The temporary directory.
_temp_dir = None


def _clear_temp_dir():
    """Clean up the temporary directory."""
    global _temp_dir
    if _temp_dir:
        tmp = _temp_dir
        _temp_dir = None
        shutil.rmtree(tmp)
    return None

atexit.register(_clear_temp_dir)


def make_temp_dir(dir=None):
    """Return the path to a temporary directory."""
    global _temp_dir
    if not _temp_dir:
        name = os.path.splitext(os.path.basename(__file__))[0]
        prefix = datetime.datetime.now().strftime(
            '{0}.tmp.%Y%m%d%H%M%S.'.format(name))
        _temp_dir = tempfile.mkdtemp(prefix=prefix, dir=dir)
    return _temp_dir


def generate_turbo(graph, out):
    """Generate turbo rules for the given graph."""
    # Skip the trivial diagram.
    if graph.ecount() == 1:
        return

    # The topology name.
    toponame = graph['name']
    # The number of denominators.
    ndenom = graph.ecount() - 2
    # The dot products as irreducible numerators.
    dots = graph['dots']
    # The non-integer parts of indices.
    nonints = graph.es['nonint'][1:-1]

    # The file name of the FORM program.
    frm_name = os.path.join(make_temp_dir(), toponame + ".frm")

    # dict for: dependent momentum -> sum of independent momenta
    rels = pathforcer.solvePropagatorSystem(graph)[2]

    # The edge IDs for independent momenta.
    st = graph.spanning_tree(return_tree=False)
    ll = [i for i in range(graph.ecount()) if i not in st]

    # The edge IDs for the external lines.
    i1 = next(e.index for e in graph.es if 0 in e.tuple)
    i2 = next(e.index for e in graph.es if graph.vcount() - 1 in e.tuple)

    # Construct parameters to be passed to FORM.
    internals = ['p{0}'.format(i) for i in ll]
    externals = ['p0']
    propagators = ['p{0}.p{0}'.format(i)
                   for i in range(graph.ecount())
                   if i != i1 and i != i2]
    propagators += ['p{0}.p{1}'.format(i, j) for i, j in dots]
    replacements = ['id {0} = {1};'.format(lhs, rhs)
                    for lhs, rhs in rels.items()]
    nonintshifts = ['id n{0} = n{0} + {1} * ep;'.format(i + 1, nonints[i])
                    for i in range(0, len(nonints))
                    if nonints[i] != 0]

    internals = ','.join(internals)
    externals = ','.join(externals)
    propagators = ','.join(propagators)
    replacements = ' '.join(replacements)
    nonintshifts = ' '.join(nonintshifts)

    with open(frm_name, 'w') as f:
        f.write('''#-
#include ibp.h
Format nospaces;
Format 120;
CF rat;
PolyRatFun rat;
V p0,...,p19;
AutoDeclare S x,y,n,i,M;
S ep,Q2,D;
CF Z;
#define Internals "{1}"
#define Externals "{2}"
#define Propagators "{3}"
#define Replacements "{4}"
#call GenerateIBP(A,Z{0},n,D)
id many, x? = rat(x,1);
multiply replace_(D,4-2*ep);
id p0.p0 = Q2;
.sort
argument rat;
  {5}
endargument;
.sort

#define N "`NumberOfPropagators'"
#define NUMEQ "`NumberOfRelations'"

#do i=1,{6}
    #$nsign`i' = 1;
    #define nsign`i' "+"
#enddo
#do i={6}+1,`N'
    #$nsign`i' = -1;
    #define nsign`i' "-"
#enddo

id  Z{0}(?x,x1?,...,x`N'?) =
    Z(<(n1-x1)*$nsign1>+...+<(n`N'-x`N')*$nsign`N'>,x1,...,x`N');
.sort
#define TOPO "{0}"
#define ORDER "`N',1,-1"

.sort
#do i = 1,`NUMEQ'
G   B`i' = A`i';
#enddo
.sort

#call NewCleanup(A,B,`N',3)
#call BuildTurbo({0},B,`N',stat.tmp)
.end
'''.format(
            toponame,
            internals,
            externals,
            propagators,
            replacements,
            nonintshifts,
            ndenom)
        )
    try:
        subprocess.check_output(['form', frm_name], universal_newlines=True)
    except subprocess.SubprocessError as e:
        subprocess.call(['cat', frm_name])
        print(e.output)
        raise

    with open('stat.tmp', 'r') as f:
        lines = f.read()
        if lines:
            out.write('''\
#do ired = 1,1
#$ired = 0;
{0}goto 2;
label 1;
id Q2 = `Q2';
$ired=1;
label 2;
ModuleOption,maximum,$ired;
.sort:turbo-{1};
#if ( `$ired\' > 0 )
\t#redefine ired "0"
#endif
#enddo
'''.format(lines, toponame))


def convert_tree_to_turbo(metagraph, destination='turbo.h'):
    """Convert a topology tree to turbo rules."""
    print('Generating "{0}"...'.format(destination))
    with open(destination, 'w') as f:
        for v in range(metagraph.vcount()):
            proc_name = 'TurboM{0}'.format(metagraph.vs[v]['graph']['name'])
            f.write('*--#[ {0} :\n'.format(proc_name))
            f.write("#procedure {0}()\n".format(proc_name))
            generate_turbo(metagraph.vs[v]['graph'], f)
            f.write('#endprocedure\n')
            f.write('*--#] {0} : \n'.format(proc_name))

    print('Done!')

if __name__ == '__main__':
    if len(sys.argv) < 2:
        print('Usage: turbo.py <nloops>')
        exit(1)

    nloops = int(sys.argv[1])

    mgfile = 'metagraph{0}.pickle'.format(nloops)

    print('Loading from "{0}"...'.format(mgfile))
    with open('metagraph{0}.pickle'.format(nloops), 'rb') as f:
        mgtree = pickle.load(f)

    convert_tree_to_turbo(mgtree)
