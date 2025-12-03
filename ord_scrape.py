import base64
import json
import os
import re
import sys
import time
from typing import Dict, List, Optional

import requests
from requests.adapters import HTTPAdapter
from urllib3.util import Retry


API_BASE = "https://open-reaction-database.org/api"
USER_AGENT = (
    "Mozilla/5.0 (Windows NT 10.0; Win64; x64) "
    "AppleWebKit/537.36 (KHTML, like Gecko) "
    "Chrome/120.0.0.0 Safari/537.36 ORD-Scraper/1.1"
)


def make_session() -> requests.Session:
    session = requests.Session()
    retries = Retry(
        total=5,
        connect=5,
        read=5,
        backoff_factor=0.8,
        status_forcelist=[429, 500, 502, 503, 504],
        allowed_methods=["GET", "POST"],
        raise_on_status=False,
    )
    adapter = HTTPAdapter(max_retries=retries)
    session.mount("http://", adapter)
    session.mount("https://", adapter)
    session.headers.update({
        "User-Agent": USER_AGENT,
        "Accept": "application/json, text/plain, */*",
    })
    return session


def get_json(session: requests.Session, url: str, params: Optional[Dict] = None) -> any:
    resp = session.get(url, params=params, timeout=30)
    resp.raise_for_status()
    return resp.json()


def submit_query(session: requests.Session, dataset_id: str, limit: int = 100) -> str:
    url = f"{API_BASE}/submit_query"
    params = {"dataset_id": dataset_id, "limit": limit}
    resp = session.get(url, params=params, timeout=30)
    resp.raise_for_status()
    task_id = resp.text.strip().strip('"')
    return task_id


def fetch_query_result(session: requests.Session, task_id: str, poll_timeout_s: int = 60) -> List[Dict]:
    url = f"{API_BASE}/fetch_query_result"
    params = {"task_id": task_id}
    start = time.time()
    while True:
        resp = session.get(url, params=params, timeout=30)
        if resp.status_code == 200:
            try:
                return resp.json()
            except json.JSONDecodeError:
                return []
        elif resp.status_code in (202, 400):
            if time.time() - start > poll_timeout_s:
                raise TimeoutError(f"Query result not ready for task_id={task_id} within {poll_timeout_s}s (status={resp.status_code}).")
            time.sleep(1.0)
            continue
        else:
            resp.raise_for_status()


def decode_reaction_proto(proto_b64: str):
    try:
        from ord_schema.proto import reaction_pb2
    except Exception:
        print("ERROR: ord_schema is not installed. Install with: pip install ord-schema", file=sys.stderr)
        raise
    raw = base64.b64decode(proto_b64)
    rxn = reaction_pb2.Reaction()
    rxn.ParseFromString(raw)
    return rxn


def compound_identifiers_to_text(identifiers) -> str:
    texts: List[str] = []
    for ident in identifiers:
        try:
            val = ident.value
            if not val:
                continue
            texts.append(val)
        except Exception:
            continue
    seen = set()
    out = []
    for x in texts:
        if x not in seen:
            out.append(x)
            seen.add(x)
    return "; ".join(out)


CORE_KEYS = {"base", "solvent", "amine", "aryl halide", "metal", "ligand"}
EXTRA_WHITELIST = {"carboxylic acid", "additive", "activation agent"}


def normalize_key(k: str) -> str:
    k = k.strip().lower().replace("_", " ")
    k = k.replace("aryl halides", "aryl halide").replace("arylhalide", "aryl halide")
    return k


def classify_target(norm: str) -> Optional[str]:
    if any(x in norm for x in ["base"]):
        return "base"
    if "solvent" in norm:
        return "solvent"
    if "amine" in norm:
        return "amine"
    if "aryl" in norm and ("halide" in norm or "bromide" in norm or "chloride" in norm or "iodide" in norm):
        return "aryl halide"
    if ("metal" in norm) or ("catalyst" in norm):
        return "metal"
    if "ligand" in norm:
        return "ligand"
    if norm in EXTRA_WHITELIST:
        return norm
    if re.fullmatch(r"m\d(?:_m\d)*", norm):
        return norm
    return None


def map_role_enum(role_enum) -> Optional[str]:
    try:
        from ord_schema.proto import reaction_pb2
        try:
            return reaction_pb2.ReactionRole.Name(role_enum)
        except Exception:
            pass
    except Exception:
        pass
    try:
        numeric = int(role_enum)
    except Exception:
        numeric = None
    if numeric is not None:
        return {1: "REACTANT", 2: "REAGENT", 3: "SOLVENT", 4: "CATALYST"}.get(numeric, str(numeric))
    return str(role_enum)


def extract_roles(rxn) -> Dict[str, Dict]:
    roles_values: Dict[str, Optional[str]] = {k: None for k in CORE_KEYS}
    roles_labels: Dict[str, Optional[str]] = {k: None for k in CORE_KEYS}
    extras_values: Dict[str, Optional[str]] = {}
    extras_labels: Dict[str, Optional[str]] = {}
    roles_raw: Dict[str, List[Dict]] = {k: [] for k in CORE_KEYS}

    try:
        inputs_map = rxn.inputs
    except Exception:
        inputs_map = {}

    for key in inputs_map:
        norm = normalize_key(key)
        target = classify_target(norm)
        if not target:
            continue
        is_core = target in CORE_KEYS
        try:
            comp_texts: List[str] = []
            comp_roles: List[str] = []
            reaction_input = inputs_map[key]
            for comp in reaction_input.components:
                text = compound_identifiers_to_text(comp.identifiers)
                if text:
                    comp_texts.append(text)
                try:
                    role_name = map_role_enum(comp.reaction_role)
                except Exception:
                    role_name = None
                if role_name:
                    comp_roles.append(role_name)
                if is_core:
                    roles_raw[target].append({"value": text, "reaction_role": role_name})
                else:
                    roles_raw.setdefault(target, []).append({"value": text, "reaction_role": role_name})
            comp_texts = [t for t in comp_texts if t]
            comp_roles = [r for r in comp_roles if r]
            if comp_texts:
                joined = "; ".join(sorted(set(comp_texts)))
                if is_core:
                    roles_values[target] = joined
                else:
                    extras_values[target] = joined
            if comp_roles:
                joinedr = "; ".join(sorted(set(comp_roles)))
                if is_core:
                    roles_labels[target] = joinedr
                else:
                    extras_labels[target] = joinedr
        except Exception:
            continue

    return {"values": roles_values, "roles": roles_labels, "extras_values": extras_values, "extras_roles": extras_labels, "raw": roles_raw}


def scrape_ord_browse(max_datasets: Optional[int] = None, per_dataset_limit: int = 0, dataset_ids: Optional[List[str]] = None) -> List[Dict]:
    session = make_session()
    datasets = get_json(session, f"{API_BASE}/datasets")
    results: List[Dict] = []
    count = 0
    dsid_set = set(dataset_ids or [])
    for ds in datasets:
        dataset_id = ds.get("dataset_id")
        name = ds.get("name")
        description = ds.get("description")
        num_rxns = ds.get("num_reactions")

        if dsid_set and dataset_id not in dsid_set:
            continue
        count += 1
        if max_datasets and count > max_datasets:
            break

        try:
            effective_limit = per_dataset_limit if per_dataset_limit and per_dataset_limit > 0 else (num_rxns or 100000)
            task_id = submit_query(session, dataset_id, limit=effective_limit)
            items = fetch_query_result(session, task_id, poll_timeout_s=60)
        except Exception as e:
            results.append({
                "dataset_id": dataset_id,
                "raw": {},
                "error": str(e),
            })
            continue

        for item in items:
            reaction_id = item.get("reaction_id")
            proto_b64 = item.get("proto")
            if not proto_b64:
                continue
            try:
                rxn = decode_reaction_proto(proto_b64)
                roles = extract_roles(rxn)
            except Exception as e:
                roles = {"values": {k: None for k in CORE_KEYS}, "roles": {k: None for k in CORE_KEYS}, "extras_values": {}, "extras_roles": {}, "raw": {}, "error": f"decode_failed: {e}"}

            results.append({
                "dataset_id": dataset_id,
                "raw": roles.get("raw", {}),
            })

            # Flatten extras into top-level fields
            try:
                extras_v = roles.get("extras_values", {})
                extras_r = roles.get("extras_roles", {})
                for k, v in extras_v.items():
                    results[-1][k] = v
                for k, r in extras_r.items():
                    results[-1][f"{k}_role"] = r
            except Exception:
                pass

        time.sleep(0.1)

    return results


def write_json(rows: List[Dict], json_path: str) -> None:
    with open(json_path, "w", encoding="utf-8") as f:
        json.dump(rows, f, ensure_ascii=False, indent=2)


def main():
    import argparse
    parser = argparse.ArgumentParser(description="Scrape ORD browse API endpoints and extract reaction component roles (JSON only).")
    parser.add_argument("--all", action="store_true")
    parser.add_argument("--max_datasets", type=int, default=None)
    parser.add_argument("--per_dataset_limit", type=int, default=0)
    parser.add_argument("--json_out", default=os.path.join(os.getcwd(), "ord_browse.json"))
    parser.add_argument("--dataset_ids", type=str, default=None)
    parser.add_argument("--reactants_only", action="store_true")
    parser.add_argument("targets", nargs="*", default=[])
    args = parser.parse_args()

    ds_ids: Optional[List[str]] = None
    if args.dataset_ids:
        ds_ids = [x.strip() for x in args.dataset_ids.split(",") if x.strip()]
    extracted_ids: List[str] = []
    for t in args.targets:
        m = re.search(r"ord_dataset-[0-9a-f]+", t)
        if m:
            extracted_ids.append(m.group(0))
        elif t.startswith("ord_dataset-"):
            extracted_ids.append(t)
    if extracted_ids:
        ds_ids = extracted_ids
    if args.all:
        args.max_datasets = None
        args.per_dataset_limit = 0

    try:
        rows = scrape_ord_browse(max_datasets=args.max_datasets, per_dataset_limit=args.per_dataset_limit, dataset_ids=ds_ids)
    except Exception as e:
        print(f"Scrape failed: {e}", file=sys.stderr)
        sys.exit(1)

    if args.reactants_only:
        filtered: List[Dict] = []
        for row in rows:
            raw = row.get("raw") or {}
            new_raw = {}
            for k, arr in raw.items():
                try:
                    new_raw[k] = [x for x in arr if (str(x.get("reaction_role")) == "REACTANT")]
                except Exception:
                    new_raw[k] = []
            filtered.append({"dataset_id": row.get("dataset_id"), "raw": new_raw})
        rows = filtered

    write_json(rows, args.json_out)
if __name__ == "__main__":
    main()
