#!/usr/bin/env python3
"""Reconstruct file states at the 'Very very very good' checkpoint by replaying
Write/Edit tool operations from the conversation transcript up to that message."""
import json, os, re

JSONL = ("/sessions/eloquent-wonderful-cray/mnt/.claude/projects/"
         "C--Users-alain-frances-AppData-Roaming-Claude-local-agent-mode-sessions-"
         "78b769ff-7fe6-4ec4-84c9-eb29c0477bf7-829fbbc7-4360-43e9-818d-ad34c10ee3f8-"
         "local-f563e60c-6806-4b9a-b5b1-ac43b5af0927-outputs/"
         "891cdce7-8ae1-4ed7-9a85-ed9d86db2215.jsonl")

CHECKPOINT_UUID = "a2c92852-5532-4656-a808-10b48f266dd4"

files = {}
read_seed = {}
tool_by_id = {}
applied = {"write": 0, "edit": 0, "edit_miss": 0, "seed": 0}


def strip_read(text):
    out, hit = [], 0
    pat = re.compile(r"^\s*\d+\t(.*)$")
    for ln in text.split("\n"):
        m = pat.match(ln)
        if m:
            out.append(m.group(1)); hit += 1
    return "\n".join(out) if hit else None


def norm_path(p):
    p = (p or "").replace("\\", "/")
    for key in ("/trunk/", "/tests/"):
        i = p.lower().find(key)
        if i != -1:
            return p[i + 1:]
    return None


with open(JSONL, encoding="utf-8") as fh:
    for line in fh:
        line = line.strip()
        if not line:
            continue
        try:
            rec = json.loads(line)
        except Exception:
            continue
        if rec.get("uuid") == CHECKPOINT_UUID:
            break
        msg = rec.get("message")
        if not isinstance(msg, dict):
            continue
        content = msg.get("content")
        if not isinstance(content, list):
            continue
        for blk in content:
            if not isinstance(blk, dict):
                continue
            t = blk.get("type")
            if t == "tool_use":
                tool_by_id[blk.get("id")] = (blk.get("name"), blk.get("input", {}))
                name, inp = blk.get("name"), blk.get("input", {})
                rel = norm_path(inp.get("file_path", "")) if inp else None
                if name == "Write" and rel:
                    files[rel] = inp.get("content", ""); applied["write"] += 1
                elif name in ("Edit", "MultiEdit") and rel:
                    edits = ([{"old_string": inp.get("old_string", ""),
                              "new_string": inp.get("new_string", ""),
                              "replace_all": inp.get("replace_all", False)}]
                             if name == "Edit" else inp.get("edits", []))
                    if rel not in files:
                        files[rel] = read_seed.get(rel, "")
                        if rel in read_seed:
                            applied["seed"] += 1
                    for e in edits:
                        old, new = e.get("old_string", ""), e.get("new_string", "")
                        if old and old in files[rel]:
                            files[rel] = (files[rel].replace(old, new)
                                          if e.get("replace_all") else
                                          files[rel].replace(old, new, 1))
                            applied["edit"] += 1
                        else:
                            applied["edit_miss"] += 1
            elif t == "tool_result":
                nm, inp = tool_by_id.get(blk.get("tool_use_id"), (None, None))
                if nm == "Read" and inp:
                    rel = norm_path(inp.get("file_path", ""))
                    if rel and rel not in read_seed:
                        c = blk.get("content")
                        if isinstance(c, list):
                            c = "\n".join(x.get("text", "") for x in c
                                          if isinstance(x, dict))
                        raw = strip_read(c) if isinstance(c, str) else None
                        if raw:
                            read_seed[rel] = raw

# --- second pass: which files were edited AFTER the checkpoint? -------------
edited_after = set()
seen_checkpoint = False
_tool = {}
with open(JSONL, encoding="utf-8") as fh:
    for line in fh:
        line = line.strip()
        if not line:
            continue
        try:
            rec = json.loads(line)
        except Exception:
            continue
        if rec.get("uuid") == CHECKPOINT_UUID:
            seen_checkpoint = True
        if not seen_checkpoint:
            continue
        msg = rec.get("message")
        if not isinstance(msg, dict):
            continue
        content = msg.get("content")
        if not isinstance(content, list):
            continue
        for blk in content:
            if isinstance(blk, dict) and blk.get("type") == "tool_use" \
                    and blk.get("name") in ("Write", "Edit", "MultiEdit"):
                rel = norm_path((blk.get("input") or {}).get("file_path", ""))
                if rel:
                    edited_after.add(rel)

BASE = "/sessions/eloquent-wonderful-cray/mnt/tmp_claude_marmites/MARMITES"
OUT = os.path.join(BASE, "trunk/MM_MF6_conversion")
n = 0
for rel, content in sorted(files.items()):
    disk = os.path.join(BASE, rel)
    # A file NOT edited after the checkpoint is unchanged on disk since then,
    # so the current disk copy IS the authoritative checkpoint version. Only
    # files edited afterwards need the (exact) replay reconstruction.
    source = "replay"
    if rel not in edited_after and os.path.exists(disk):
        with open(disk, encoding="utf-8", errors="replace") as f:
            content = f.read()
        source = "disk "
    dest = os.path.join(OUT, rel)
    os.makedirs(os.path.dirname(dest), exist_ok=True)
    with open(dest, "w", encoding="utf-8") as f:
        f.write(content)
    n += 1
    print("  [%s] %-50s %6d bytes" % (source, rel, len(content)))
print("\nreconstructed %d files" % n)
print("ops:", applied)
