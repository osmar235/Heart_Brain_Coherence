import os, sys, re, json, hashlib, csv, collections, statistics as st
repo, annex, out = sys.argv[1], sys.argv[2], sys.argv[3]
UUID_ON = "8e7ce7f9-9140-4bd4-99bd-d550bca0e4de"
# key -> (versionId, s3path), presence
ver = {}; present = {}
for root, d, files in os.walk(annex):
    if ".git" in root.split(os.sep): continue
    for f in files:
        p = os.path.join(root, f)
        if f.endswith(".log.rmet"):
            key = f[:-9]
            for line in open(p, encoding="utf-8"):
                m = re.match(r"\S+ (\S+):V \+([^#]+)#(.+)$", line.strip())
                if m and m.group(1) == UUID_ON: ver.setdefault(key, {})[m.group(3)] = m.group(2)
        elif f.endswith(".log") and f.startswith("SHA256E"):
            key = f[:-4]
            for line in open(p, encoding="utf-8"):
                a = line.split()
                if len(a) == 3 and a[2] == UUID_ON: present[key] = a[1]
rows = []
for root, d, files in os.walk(repo):
    d[:] = [x for x in d if x != ".git"]
    for f in sorted(files):
        p = os.path.join(root, f); rel = os.path.relpath(p, repo).replace(os.sep, "/")
        if rel.startswith((".nemar", ".datalad", ".github")) or rel == ".gitattributes": continue
        m = re.match(r"sub-\d+", rel); sub = m.group(0) if m else "root"
        t = re.search(r"task-([a-z]+)", rel); task = t.group(1) if t else "na"
        if os.path.islink(p):
            key = os.path.basename(os.readlink(p))
            k = re.match(r"SHA256E-s(\d+)--([0-9a-f]{64})", key)
            vv = ver.get(key, {}); vid = vv.get("ds008302/" + rel, "")
            rows.append([rel, sub, task, "openneuro_s3", k.group(1), k.group(2), vid, present.get(key, "?"), "1" if vid else "0"])
        else:
            b = open(p, "rb").read()
            rows.append([rel, sub, task, "github_raw", str(len(b)), hashlib.sha256(b).hexdigest(), "", "", ""])
rows.sort()
with open(os.path.join(out, "manifest.tsv"), "w", newline="\n", encoding="utf-8") as fh:
    w = csv.writer(fh, delimiter="\t", lineterminator="\n")
    w.writerow(["path", "subject", "task", "source", "size_bytes", "sha256", "s3_version_id", "s3_present_flag", "s3_path_matches"])
    w.writerows(rows)
print("files", len(rows), "annex", sum(r[3] == "openneuro_s3" for r in rows), "bytes", sum(int(r[4]) for r in rows))
print("annex without versionId", sum(r[3] == "openneuro_s3" and not r[6] for r in rows), "present!=1", sum(r[3] == "openneuro_s3" and r[7] != "1" for r in rows), "path mismatch", sum(r[3] == "openneuro_s3" and r[8] != "1" for r in rows))
size = {r[0]: int(r[4]) for r in rows}
# audit
subs = sorted(x for x in os.listdir(repo) if x.startswith("sub-"))
A = []
for s in subs:
    e = os.path.join(repo, s, "eeg"); r = {"sub": s}
    for t in ["memory", "rest", "video"]:
        setp = f"{s}/eeg/{s}_task-{t}_eeg.set"; r[t + "_gb"] = round(size.get(setp, 0) / 1e9, 2) if setp in size else None
        jp = os.path.join(e, f"{s}_task-{t}_eeg.json")
        if os.path.exists(jp):
            j = json.load(open(jp)); r[t + "_dur"] = j.get("RecordingDuration"); r[t + "_fs"] = j.get("SamplingFrequency"); r[t + "_ref"] = j.get("EEGReference")
        cp = os.path.join(e, f"{s}_task-{t}_channels.tsv")
        if os.path.exists(cp):
            ch = list(csv.DictReader(open(cp), delimiter="\t"))
            r[t + "_nch"] = len(ch); r[t + "_ecg"] = ";".join(f"{c['name']}/{c['type']}/{c.get('units')}" for c in ch if "ECG" in (c["name"] + c["type"]).upper()); r[t + "_ppg"] = ";".join(c["name"] for c in ch if "PPG" in c["name"].upper())
    r["beh_files"] = len(os.listdir(os.path.join(repo, s, "beh"))) if os.path.isdir(os.path.join(repo, s, "beh")) else 0
    ep = os.path.join(e, f"{s}_task-memory_events.tsv")
    if os.path.exists(ep):
        ev = [(float(x["onset"]), x["duration"], x["value"]) for x in csv.DictReader(open(ep), delimiter="\t")]
        dur = r.get("memory_dur") or ev[-1][0]
        r["interior_boundaries"] = sum(1 for o, d_, v in ev if v == "boundary" and 1 < o < dur - 1)
        r["buf_overflow"] = sum(1 for o, d_, v in ev if v == "Buffer Overflow")
        # trials
        trials = []; cur = None; blocks = []; bstart = None
        for o, d_, v in ev:
            if v == "S 90":
                if not blocks or o - blocks[-1]["t0"] > 1: blocks.append({"t0": o, "trials": [], "t1": None})
            m = re.match(r"S1(0\d|1[0-3])$", v)
            if m: cur = {"len": int(m.group(1)), "t0": o, "flags": [], "probe": None}
            dm = re.match(r"(\d\d)(\d\d)([01x])$", v)
            if dm and cur is not None: cur["flags"].append(dm.group(3))
            if v == "S160" and cur is not None:
                cur["probe"] = o
                if blocks: blocks[-1]["trials"].append(cur)
                trials.append(cur); cur = None
            if v in ("S180", "S190") and blocks and blocks[-1]["t1"] is None: blocks[-1]["t1"] = o
        ib = [(o, d_) for o, d_, v in ev if v == "boundary" and 1 < o < dur - 1]
        r["bnd_in_block"] = sum(1 for o, d_ in ib if any(b["t1"] and b["t0"] <= o <= b["t1"] for b in blocks))
        r["bnd_list"] = [[round(o,1), d_] for o, d_ in ib]
        segs = [0.0] + [o for o, d_ in ib] + [dur]
        r["longest_seg"] = round(max(segs[i+1]-segs[i] for i in range(len(segs)-1)))
        r["n_trials"] = len(trials); r["n_blocks"] = len(blocks)
        r["unknown_flags"] = sum(1 for t in trials if "x" in t["flags"])
        r["bad_len_trials"] = sum(1 for t in trials if len(t["flags"]) != t["len"])
        cnt = collections.Counter(t["len"] for t in trials); r["len_counts"] = dict(sorted(cnt.items()))
        big = sorted(L for L, c in cnt.items() if c >= 60)
        pair = [(a, b) for a in big for b in big if b == a + 1]
        if pair:
            sp, ov = pair[-1] if len(pair) == 1 else max(pair, key=lambda ab: cnt[ab[0]] + cnt[ab[1]])
            r["span"] = sp; r["n_pairs"] = len(pair)
            for nm, L in (("span", sp), ("over", ov)):
                tt = [t for t in trials if t["len"] == L and "x" not in t["flags"] and len(t["flags"]) == L]
                r[nm + "_n"] = sum(1 for t in trials if t["len"] == L)
                r[nm + "_acc"] = round(100 * sum(all(f == "1" for f in t["flags"]) for t in tt) / len(tt), 1) if tt else None
        bd = [b["t1"] - b["t0"] for b in blocks if b["t1"] and len(b["trials"]) == 10]
        r["blk_n10"] = len(bd)
        if bd: r["blk_dur_min"] = round(min(bd)); r["blk_dur_med"] = round(st.median(bd)); r["blk_dur_max"] = round(max(bd))
        br = [blocks[i + 1]["t0"] - blocks[i]["t1"] for i in range(len(blocks) - 1) if blocks[i]["t1"]]
        if br: r["break_med"] = round(st.median(br)); r["break_max"] = round(max(br))
    A.append(r)
json.dump(A, open(os.path.join(out, "audit.json"), "w"), indent=1)
