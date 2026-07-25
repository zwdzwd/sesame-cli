"""Fail a build whose linux-64 package would not run on ordinary machines.

Usage: python tools/check_glibc_floor.py <conda-build output folder>

CI runners carry a far newer glibc than users do. conda-forge's gcc_linux-64
depends on sysroot_linux-64 with no version constraint, so without an explicit
pin the solver takes the newest sysroot and the binary silently acquires a
glibc floor no ordinary machine satisfies -- e.g. `-D_GNU_SOURCE` plus glibc
>= 2.38 rewrites strtol to __isoc23_strtol@GLIBC_2.38, and glibc 2.38 gave
fmod a new symbol version too. Packages built that way install cleanly and
then die at exec time with "version `GLIBC_2.38' not found".

Two things are checked, because either alone leaves a hole:
  * every ELF in bin/ must require nothing above GLIBC_2.17 (the pin worked);
  * the package must declare a `__glibc >=` run constraint (from
    {{ stdlib('c') }}), so conda refuses to install it on a host too old to
    run it rather than failing later at exec time.

Kept in Python rather than shell: the equivalent shell version interacted
badly with GitHub's `bash -el` and failed builds that were actually fine.
"""
import glob, json, os, re, subprocess, sys, tarfile, tempfile, zipfile

FLOOR = (2, 17)
pkgs = sorted(glob.glob(os.path.join(sys.argv[1], "linux-64", "*.conda")))
if not pkgs:
    sys.exit("no .conda package found in %s/linux-64" % sys.argv[1])
pkg = pkgs[0]
print("checking", os.path.basename(pkg))

work = tempfile.mkdtemp()
with zipfile.ZipFile(pkg) as z:
    z.extractall(work)
for member in glob.glob(os.path.join(work, "*.tar.zst")):
    # python's tarfile has no zstd support before 3.14; GNU tar auto-detects it.
    for cmd in (["tar", "-xf", member, "-C", work],
                ["tar", "--use-compress-program=unzstd", "-xf", member, "-C", work]):
        if subprocess.run(cmd, capture_output=True).returncode == 0:
            break
    else:
        sys.exit("could not extract %s" % member)

idx = json.load(open(os.path.join(work, "info", "index.json")))
depends = idx.get("depends", [])
print("depends:", depends)
failures = []

# The pin must show up as a run constraint, otherwise conda will happily
# install this on a host too old to run it -- the failure mode that let the
# broken build reach users in the first place.
if not any(d.startswith("__glibc >=") for d in depends):
    failures.append("package has no __glibc run constraint; is {{ stdlib('c') }} in meta.yaml?")

bindir = os.path.join(work, "bin")
binaries = sorted(glob.glob(os.path.join(bindir, "*"))) if os.path.isdir(bindir) else []
if not binaries:
    failures.append("no files found in the package's bin/")

for b in binaries:
    with open(b, "rb") as fh:
        if fh.read(4) != b"\x7fELF":
            print("  %s: not an ELF file, skipped" % os.path.basename(b))
            continue
    out = subprocess.run(["objdump", "-T", b], capture_output=True, text=True)
    if out.returncode != 0:
        failures.append("objdump failed on %s: %s" % (b, out.stderr.strip()[:200]))
        continue
    vers = set()
    for m in re.finditer(r"GLIBC_(\d+)\.(\d+)(?:\.(\d+))?", out.stdout):
        vers.add((int(m.group(1)), int(m.group(2))))
    if not vers:
        print("  %s: references no versioned glibc symbols" % os.path.basename(b))
        continue
    top = max(vers)
    label = "%d.%d" % top
    if top > FLOOR:
        offenders = sorted({
            ln.split()[-1] for ln in out.stdout.splitlines()
            if re.search(r"GLIBC_%d\.%d\b" % top, ln)
        })
        failures.append(
            "%s requires GLIBC_%s, above the pinned %d.%d floor; symbols: %s"
            % (os.path.basename(b), label, FLOOR[0], FLOOR[1], ", ".join(offenders[:8]))
        )
    else:
        print("  %s: highest requirement GLIBC_%s -- ok" % (os.path.basename(b), label))

for f in failures:
    print("::error::" + f)
    sys.stderr.write("FAIL: %s\n" % f)
sys.exit(1 if failures else 0)
