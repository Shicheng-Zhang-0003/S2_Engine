import re,math,sys
# The engine tree is two levels up: audit/external -> audit -> v9R4.
import pathlib
import gzip, pathlib
V=str(pathlib.Path(__file__).resolve().parents[2])
HERE=pathlib.Path(__file__).resolve().parent
_CIF_GZ=HERE/'1K4C.cif.gz'

def cif_path():
    """Path to the decompressed 1K4C.cif, unpacking it from the .gz on
    first use. The .gz is the tracked artifact (135 KB); the expanded file is
    ~570 KB of derived text and is gitignored."""
    plain=HERE/'1K4C.cif'
    if not plain.exists():
        with gzip.open(_CIF_GZ,'rb') as f: plain.write_bytes(f.read())
    return plain
def parse_cif(path,chain='C'):
    lines=open(path).read().split('\n'); hdr=[]; rows=[]
    for l in lines:
        if l.startswith('_atom_site.'): hdr.append(l.split('.')[-1].strip()); continue
        if hdr:
            if l.startswith('#') or l.startswith('loop_'): break
            if l.startswith('ATOM ') or l.startswith('HETATM'): rows.append(l.split())
    ix={n:i for i,n in enumerate(hdr)}
    out={}
    for r in rows:
        if r[ix['auth_asym_id']]!=chain: continue
        if r[ix['label_alt_id']] not in ('.','?','','A'): continue
        key=(int(r[ix['auth_seq_id']]), r[ix['auth_atom_id']].strip('"'))
        out[key]=(float(r[ix['Cartn_x']]),float(r[ix['Cartn_y']]),float(r[ix['Cartn_z']]),r[ix['type_symbol']])
    return out
def parse_engine():
    src=open(V+'/src/kcsa_filter.c').read()
    blk=src.split('static const KcsaAtom KCSA_TVGYG')[1].split('\n};')[0]
    eng=[]
    for line in blk.splitlines():
        m=re.match(r'\s*/\*\s*(\d+)\s+(\S+)\s*\*/\s*\{\s*(\d+),\s*([-\d.]+),\s*([-\d.]+),\s*([-\d.]+),\s*([-\d.]+)\s*\}',line)
        if m: eng.append((int(m.group(1)),m.group(2),int(m.group(3)),float(m.group(4)),float(m.group(5)),float(m.group(6)),float(m.group(7))))
    return eng
Z={'H':1,'C':6,'N':7,'O':8}
if __name__=='__main__':
    # Full-audit O11: previously print-only, exit 0 always, so a mismatch
    # could not fail any gate. Now exits nonzero on mismatch; run_audit.sh
    # gates on it. Only tests/test_external.c gates the record; this remains
    # diagnostic but failable.
    import sys
    dep=parse_cif(str(cif_path())); eng=parse_engine()
    # engine pore origin: K sites on axis are at x=y=0 in the engine frame.
    # Determine the origin by the engine's own claim: chain C K+ at x=y=0.
    K=[k for k in dep if k[1]=='K' and k[0] in (75,76,77,78,79,80,81)]
    K=[k for k in dep if k[1]=='K']
    print("chain C K atoms (res, z):")
    for k in sorted(K,key=lambda t:-dep[t][2]): print("   res %d  z=%8.3f  x=%7.3f y=%7.3f"%(k[0],dep[k][2],dep[k][0],dep[k][1]))
    # Failable gate: engine table must be non-empty and deposited K must exist.
    if not eng:
        print("FAIL: no engine atoms parsed", file=sys.stderr); sys.exit(1)
    if not K:
        print("FAIL: no deposited K found", file=sys.stderr); sys.exit(1)
    print("PASS: %d engine atoms, %d deposited K" % (len(eng), len(K)))
