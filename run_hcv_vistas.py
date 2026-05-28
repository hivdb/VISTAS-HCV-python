#!/usr/bin/env python3
"""Generate the HCV GenBank Step 2 descriptive CSV from a FASTA file.

This is a Python port of the HCV Step 1 browser workflow:

1. parse FASTA records and make request headers unique
2. call https://hivdb.stanford.edu/hcv/graphql for sequenceAnalysis
3. build the alignment report used by generateEditableTable()
4. merge multiple gene records for the same isolate
5. run the HCV genotyping pass
6. call EASL-extend drug-resistance positions
7. write the Step 2 descriptive CSV
"""

from __future__ import annotations

import argparse
import csv
import json
import re
import sys
import time
import uuid
from collections import Counter, OrderedDict, defaultdict
from pathlib import Path
from urllib.error import HTTPError, URLError
from urllib.request import Request, urlopen


SCRIPT_DIR = Path(__file__).resolve().parent
DEFAULT_ENDPOINT = "https://hivdb.stanford.edu/hcv/graphql"

GENE_NAMES = [
    "CORE",
    "E1",
    "E2",
    "P7",
    "NS2",
    "NS3",
    "NS4A",
    "NS4B",
    "NS5A_NTD",
    "NS5A_CTD",
    "NS5B",
]
GENOTYPING_GENE_NAMES = ["NS3", "NS5A_NTD", "NS5B"]
SORT_GENE_NAMES = [
    "polyprotein",
    "CORE",
    "E1",
    "E2",
    "P7",
    "NS2",
    "NS3",
    "NS4A",
    "NS4B",
    "NS5A",
    "NS5A_NTD",
    "NS5A_CTD",
    "NS5B",
]
DRM_GENES = ["NS3", "NS5A_NTD", "NS5B"]

ALIGN_QUERY = """
query align($sequences: [UnalignedSequenceInput]!, $geneNames: [EnumGene]) {
  viewer {
    currentVersion {
      text
      publishDate
    }
    sequenceAnalysis(sequences: $sequences) {
      inputSequence {
        header
        SHA512
        sequence
      }
      bestMatchingSubtype {
        displayWithoutDistance
      }
      strain {
        name
        display
      }
      validationResults {
        level
        message
      }
      DRMs: mutations(filterOptions: [DRM]) {
        gene {
          name
          drugClasses { name }
        }
        shortText
      }
      alignedGeneSequences(includeGenes: $geneNames) {
        gene {
          name
          length
        }
        firstAA
        lastAA
        firstNA
        lastNA
        mutations {
          position
          AAs
          isInsertion
          isDeletion
          insertedNAs
          isUnusual
          isSDRM
          isDRM
          primaryType
          hasStop
          shortText
          isApobecMutation
          isAmbiguous
        }
        prettyPairwise {
          positionLine
          refAALine
          alignedNAsLine
          mutationLine
        }
      }
    }
  }
}
"""

GENOTYPING_QUERY = """
query align($sequences: [UnalignedSequenceInput]!) {
  viewer {
    currentVersion {
      text
      publishDate
    }
    sequenceAnalysis(sequences: $sequences) {
      inputSequence {
        header
        SHA512
        sequence
      }
      bestMatchingSubtype {
        displayWithoutDistance
      }
      strain {
        name
        display
      }
    }
  }
}
"""


def sort_gene_key(name: str) -> int:
    try:
        return SORT_GENE_NAMES.index(name)
    except ValueError:
        return len(SORT_GENE_NAMES)


def load_json(filename: str):
    with (SCRIPT_DIR / filename).open("r", encoding="utf-8") as handle:
        return json.load(handle)


def parse_fasta(path: Path) -> list[dict[str, str]]:
    records: list[dict[str, str]] = []
    header: str | None = None
    sequence_parts: list[str] = []

    with path.open("r", encoding="utf-8") as handle:
        for line_number, raw_line in enumerate(handle, start=1):
            line = raw_line.strip()
            if not line or line.startswith("#"):
                continue
            if line.startswith(">"):
                if header is not None:
                    records.append({"header": unique_header(header), "sequence": "".join(sequence_parts)})
                header = line[1:].strip()
                sequence_parts = []
                continue
            if header is None:
                raise ValueError(f"{path}: line {line_number} contains sequence data before a FASTA header")
            sequence_parts.append(line)

    if header is not None:
        records.append({"header": unique_header(header), "sequence": "".join(sequence_parts)})
    if not records:
        raise ValueError(f"{path}: no FASTA records found")
    return records


def unique_header(header: str) -> str:
    # Step 1 appends a random token before parsing so duplicate headers survive
    # GraphQL processing and can be merged later by the first whitespace token.
    return f"{header} {int(time.time() * 1000)}-{uuid.uuid4().hex[:12]}"


def normalized_isolate_id(isolate_id: str) -> str:
    return isolate_id.split(" ")[0]


def graphql_request(endpoint: str, query: str, variables: dict) -> dict:
    payload = json.dumps(
        {"operationName": "align", "query": query, "variables": variables}
    ).encode("utf-8")
    request = Request(
        endpoint,
        data=payload,
        headers={"Content-Type": "application/json", "Accept": "application/json"},
        method="POST",
    )
    try:
        with urlopen(request, timeout=120) as response:
            body = response.read().decode("utf-8")
    except HTTPError as exc:
        detail = exc.read().decode("utf-8", errors="replace")
        raise RuntimeError(f"GraphQL HTTP {exc.code}: {detail}") from exc
    except URLError as exc:
        raise RuntimeError(f"GraphQL request failed: {exc.reason}") from exc

    result = json.loads(body)
    if result.get("errors"):
        raise RuntimeError(f"GraphQL errors: {json.dumps(result['errors'], indent=2)}")
    if "data" not in result:
        raise RuntimeError(f"GraphQL response missing data: {body}")
    return result["data"]


def progress(message: str) -> None:
    print(message, file=sys.stderr, flush=True)


def get_gaps(deletions: list[int]) -> list[list[int]]:
    gaps: list[list[int]] = []
    group: list[int] = []
    previous: int | None = None
    for current in sorted(deletions):
        if previous is not None and current - previous != 1:
            if group:
                gaps.append([group[0], group[-1], group[-1] - group[0] + 1])
            group = []
        group.append(current)
        previous = current
    if group:
        gaps.append([group[0], group[-1], group[-1] - group[0] + 1])
    return gaps


def get_na_gaps(aa_gaps: list[list[int]]) -> list[list[int]]:
    return [[start * 3 - 2, stop * 3, length * 3] for start, stop, length in aa_gaps]


def get_na_seq(gene: dict, original_seq: str) -> dict:
    gene_aa_length = gene["gene"]["length"]
    gene_na_length = gene_aa_length * 3
    first_na = gene["firstNA"]
    last_na = gene["lastNA"]
    gene_cut = original_seq[first_na - 1:last_na]

    aa_positions = [
        int(pos)
        for pos in gene["prettyPairwise"]["positionLine"]
        if str(pos).strip() != ""
    ]
    trimmed_seq_na = list(gene["prettyPairwise"]["alignedNAsLine"])

    if trimmed_seq_na and "-" in trimmed_seq_na[0]:
        trimmed_seq_na = trimmed_seq_na[1:]
        if aa_positions:
            aa_positions = aa_positions[1:]
    if trimmed_seq_na and "-" in trimmed_seq_na[-1]:
        trimmed_seq_na = trimmed_seq_na[:-1]
        if aa_positions:
            aa_positions = aa_positions[:-1]

    joined_trimmed_na = "".join(trimmed_seq_na)
    if aa_positions:
        start_aa_pos = min(aa_positions)
        stop_aa_pos = max(aa_positions)
    else:
        start_aa_pos = gene.get("firstAA") or 1
        stop_aa_pos = gene.get("lastAA") or start_aa_pos

    start_na_pos = start_aa_pos * 3 - 2
    stop_na_pos = len(joined_trimmed_na) + start_na_pos - 1
    aligned_na_length = len(joined_trimmed_na)
    aligned_na = "." * (start_na_pos - 1) + joined_trimmed_na + "." * max(gene_na_length - stop_na_pos, 0)

    untrans_reason = ""
    if start_na_pos % 3 != 1:
        untrans_reason = "start NA pos not aligned to codon"
    elif aligned_na_length % 3 != 0:
        untrans_reason = "stop NA pos not aligned to codon"
    elif len(joined_trimmed_na) % 3 != 0:
        untrans_reason = "untranslatable"

    leading = len(aligned_na) - len(aligned_na.lstrip("."))
    trailing = len(aligned_na) - len(aligned_na.rstrip("."))

    return {
        "gene_AA_length": gene_aa_length,
        "gene_NA_length": gene_na_length,
        "start_AA_pos": start_aa_pos,
        "start_NA_pos": start_na_pos,
        "stop_AA_pos": stop_aa_pos,
        "stop_NA_pos": stop_na_pos,
        "gene_AA_range": f"{start_aa_pos} ~ {stop_aa_pos}",
        "N in seq": "yes" if "." in joined_trimmed_na else "",
        "del in seq": "yes" if "-" in joined_trimmed_na else "",
        "translatable": 0 if untrans_reason else 1,
        "untrans_reason": untrans_reason,
        "trimmed_NA": joined_trimmed_na,
        "aligned_NA": aligned_na,
        "aligned_NA_heading": leading,
        "aligned_NA_heading_pos": leading + 1,
        "aligned_NA_tailing": trailing,
        "aligned_NA_tailing_pos": len(aligned_na.rstrip(".")),
        "geneCut": gene_cut,
        "firstNA": first_na,
        "lastNA": last_na,
    }


def get_genes_drms(seq: dict) -> dict[str, dict[str, list[str]]]:
    gene_classes = {"NS3": ["PI"], "NS5A_NTD": ["NS5A_I"], "NS5B": ["NI", "NNI"]}
    drms_by_gene: dict[str, list[dict]] = defaultdict(list)
    for drm in seq.get("DRMs") or []:
        gene_name = drm.get("gene", {}).get("name")
        if gene_name:
            drms_by_gene[gene_name].append(drm)

    result = {}
    for gene_name, drug_classes in gene_classes.items():
        result[gene_name] = {}
        for drug_class in drug_classes:
            mutations = []
            for drm in drms_by_gene.get(gene_name, []):
                drm_classes = [d["name"] for d in drm.get("gene", {}).get("drugClasses", [])]
                if drug_class in drm_classes:
                    mutations.append(drm.get("shortText", ""))
            result[gene_name][drug_class] = [m for m in mutations if m]
    return result


def build_alignment_report(data: dict) -> list[dict]:
    report = []
    for seq in data["viewer"]["sequenceAnalysis"]:
        original_seq = seq["inputSequence"]["sequence"]
        genes_drms = get_genes_drms(seq)
        records = []

        for gene in seq.get("alignedGeneSequences") or []:
            gene_name = gene["gene"]["name"]
            gene_mutations = gene.get("mutations") or []
            mutations = [mut["shortText"] for mut in gene_mutations if mut.get("shortText")]
            unusuals = [mut["shortText"] for mut in gene_mutations if mut.get("isUnusual") and mut.get("shortText")]
            apobecs = [mut["shortText"] for mut in gene_mutations if mut.get("isApobecMutation") and mut.get("shortText")]
            ambiguous = [mut["shortText"] for mut in gene_mutations if mut.get("isAmbiguous") and mut.get("shortText")]
            insertions = [mut["shortText"] for mut in gene_mutations if mut.get("isInsertion") and mut.get("shortText")]
            deletions = [mut["position"] for mut in gene_mutations if mut.get("isDeletion")]
            stops = [mut["shortText"] for mut in gene_mutations if mut.get("hasStop") and mut.get("shortText")]
            gaps = get_gaps(deletions)

            class_drms = genes_drms.get(gene_name, {})
            drm_text = ", ".join(
                mutation
                for drug_class in class_drms.values()
                for mutation in drug_class
            )

            alignment = get_na_seq(gene, original_seq)
            record = {
                "geneName": gene_name,
                "geneOrder": GENE_NAMES.index(gene_name) if gene_name in GENE_NAMES else len(GENE_NAMES),
                "mutations": ", ".join(mutations),
                "num_mutations": len(mutations),
                "DRMs": drm_text,
                "unusuals": ", ".join(unusuals),
                "num_unusual": len(unusuals),
                "apobecs": ", ".join(apobecs),
                "num_apobecs": len(apobecs) if apobecs else "",
                "num_ambiguous": len(ambiguous),
                "apobec_unusual": " / ".join([", ".join(apobecs), ", ".join(unusuals)]),
                "insertions": ", ".join(insertions),
                "num_insertions": len(insertions),
                "deletions": ", ".join(str(d) for d in deletions),
                "num_deletions": len(deletions),
                "gaps": "\n".join(f"[{start}-{stop}] ({length})" for start, stop, length in gaps),
                "NAgaps_tailing": [
                    gene["lastAA"] * 3 + 1 if gene["lastAA"] < gene["gene"]["length"] else gene["gene"]["length"] * 3,
                    gene["gene"]["length"] * 3,
                    (gene["gene"]["length"] - gene["lastAA"]) * 3,
                ],
                "NAgaps_tailing_exists": gene["gene"]["length"] > gene["lastAA"],
                "AAgaps": gaps,
                "NAgaps": get_na_gaps(gaps),
                "stops": ", ".join(stops),
                "alignment": alignment,
            }
            if any(gap[-1] % 3 != 0 for gap in record["NAgaps"]):
                record["alignment"]["translatable"] = False
                record["alignment"]["untrans_reason"] = "untranslatable"
            records.append(record)

        gene_info = {record["geneName"]: record for record in records}
        isolate = {
            "IsolateID": seq["inputSequence"]["header"],
            "originalSeq": original_seq,
            "virus": seq.get("strain", {}).get("name", ""),
            "subtype": (seq.get("bestMatchingSubtype") or {}).get("displayWithoutDistance", "") or "",
            "strain": seq.get("strain", {}).get("display", "") or "",
            "genes": gene_info,
            "sequenceAnalysis": seq,
        }
        report.append(isolate)

    return cut_short_gene(report)


def cut_short_gene(isolates: list[dict]) -> list[dict]:
    filtered = []
    for isolate in isolates:
        kept = {}
        for gene_name, gene in isolate["genes"].items():
            if gene_name == "NS3":
                if len(gene["alignment"]["trimmed_NA"]) > 200:
                    kept[gene_name] = gene
                continue

            if gene["num_deletions"] > gene["alignment"]["gene_AA_length"] / 2:
                continue
            if len(gene["alignment"]["geneCut"]) < gene["alignment"]["gene_AA_length"] / 2:
                continue
            if len(gene["alignment"]["trimmed_NA"]) > 150:
                kept[gene_name] = gene
        isolate["genes"] = kept
        if kept:
            filtered.append(isolate)
    return filtered


def validate_gene_exist(alignments: list[dict]) -> None:
    invalid = [iso["IsolateID"] for iso in alignments if not any(g in iso["genes"] for g in DRM_GENES)]
    if invalid:
        raise ValueError("no DRM gene found for isolate(s): " + ", ".join(invalid))


def validate_no_conflict(alignments: list[dict]) -> None:
    conflicts = [iso["IsolateID"] for iso in alignments if "conflict" in iso["IsolateID"]]
    if conflicts:
        raise ValueError("duplicated IsolateID after parsing: " + ", ".join(conflicts))


def most_common_value(values: list[str]) -> str:
    values = [value.strip() for value in values if value and value.strip()]
    preferred = [value for value in values if value.lower() != "unknown"]
    values = preferred or values
    if not values:
        return ""
    return Counter(values).most_common(1)[0][0]


def merge_isolates(isolates: list[dict]) -> list[dict]:
    grouped: OrderedDict[str, list[dict]] = OrderedDict()
    for iso in isolates:
        grouped.setdefault(normalized_isolate_id(iso["IsolateID"]), []).append(iso)

    merged = []
    for isolate_id, iso_list in grouped.items():
        main = iso_list[0]
        main["IsolateID"] = isolate_id
        if len(iso_list) > 1:
            for iso in iso_list:
                main["genes"].update(iso["genes"])
            main["subtype"] = most_common_value([iso.get("subtype", "") for iso in iso_list])
            main["strain"] = most_common_value([iso.get("strain", "") for iso in iso_list])

            seen_genes = set()
            merged_gene_sequences = []
            for iso in iso_list:
                for gene_seq in iso.get("sequenceAnalysis", {}).get("alignedGeneSequences", []) or []:
                    gene_name = gene_seq.get("gene", {}).get("name")
                    if gene_name and gene_name not in seen_genes:
                        seen_genes.add(gene_name)
                        merged_gene_sequences.append(gene_seq)
            main["sequenceAnalysis"]["alignedGeneSequences"] = merged_gene_sequences
            main.pop("originalSeq", None)
        merged.append(main)
    return merged


def run_genotyping(batch: list[dict], endpoint: str) -> list[dict]:
    partial_gene_records = []
    for seq_report in batch:
        parts = []
        for gene_name in ["NS3", "NS4A", "NS4B", "NS5A_NTD", "NS5A_CTD", "NS5B"]:
            gene = seq_report["genes"].get(gene_name)
            if not gene:
                parts.append("." * 100)
            elif gene_name in GENOTYPING_GENE_NAMES:
                parts.append(gene["alignment"]["geneCut"])
            else:
                parts.append("." * len(gene["alignment"]["geneCut"]))
        partial_gene_records.append({"header": seq_report["IsolateID"], "sequence": "".join(parts)})

    data = graphql_request(endpoint, GENOTYPING_QUERY, {"sequences": partial_gene_records})
    genotyping = {}
    for seq in data["viewer"]["sequenceAnalysis"]:
        genotyping[seq["inputSequence"]["header"]] = {
            "subtype": (seq.get("bestMatchingSubtype") or {}).get("displayWithoutDistance", "") or "",
            "strain": seq.get("strain", {}).get("display", "") or "",
        }

    for seq_report in batch:
        next_values = genotyping.get(seq_report["IsolateID"])
        if next_values:
            subtype = next_values["subtype"]
            strain = next_values["strain"]
            if subtype and subtype.strip().lower() != "unknown":
                seq_report["subtype"] = subtype
            if strain and strain.strip().lower() != "unknown":
                seq_report["strain"] = strain
        subtype = seq_report.get("subtype") or ""
        genotype_match = re.match(r"^Genotype\s+\d+", subtype, re.I)
        seq_report["genotype"] = genotype_match.group(0) if genotype_match else seq_report.get("strain", "").replace("Hepatitis C", "").strip()
    return batch


def get_genotype_from_subtype(subtype: str, fallback_display: str = "") -> str | None:
    for source in (subtype or "", fallback_display or ""):
        match = re.search(r"Genotype\s+(\d+[a-zA-Z]*)", source, re.I)
        if match:
            return match.group(1)
    return None


def get_reference_sequence(refs: list[dict], subtype: str, gene_name: str, fallback_display: str = "") -> str | None:
    genotype = get_genotype_from_subtype(subtype, fallback_display)
    if not genotype:
        return None
    normalized = genotype.lower()
    if normalized in {"1a", "1b"}:
        strain = f"hcv{normalized}"
        for ref in refs:
            if ref.get("strain", "").lower() == strain and ref.get("gene") == gene_name:
                return ref.get("refSequence")
    genotype_number = re.match(r"^(\d+)", genotype)
    if not genotype_number:
        return None
    strain = f"HCV{genotype_number.group(1)}"
    for ref in refs:
        if ref.get("strain") == strain and ref.get("gene") == gene_name:
            return ref.get("refSequence")
    return None


def grouped_pretty_pairwise_mutations(gene: dict) -> list[dict]:
    pp = gene.get("prettyPairwise") or {}
    position_line = pp.get("positionLine") or []
    ref_aa_line = pp.get("refAALine") or []
    mutation_line = pp.get("mutationLine") or []
    grouped = []
    current = None
    last_position = None

    for idx, raw_position in enumerate(position_line):
        raw = str(raw_position or "").strip()
        parsed_position = int(raw) if raw else last_position
        if parsed_position is None:
            continue
        last_position = parsed_position
        if current is None or current["position"] != parsed_position:
            if current is not None:
                grouped.append(current)
            current = {"position": parsed_position, "mutation_parts": [], "ref_aa_parts": []}
        mutation_aa = str(mutation_line[idx] if idx < len(mutation_line) else "").strip()
        ref_aa = str(ref_aa_line[idx] if idx < len(ref_aa_line) else "").strip()
        if mutation_aa:
            current["mutation_parts"].append(mutation_aa)
        if ref_aa:
            current["ref_aa_parts"].append(ref_aa)
    if current is not None:
        grouped.append(current)

    result = []
    for item in grouped:
        ref_joined = "".join(item["ref_aa_parts"])
        fallback_match = re.search(r"[A-Za-z*]", ref_joined)
        result.append(
            {
                "position": item["position"],
                "mutAA": "".join(item["mutation_parts"]),
                "fallbackRefAA": fallback_match.group(0) if fallback_match else "",
            }
        )
    return result


def normalize_genotype(genotype: str) -> str:
    if genotype in {"1a", "1b"}:
        return genotype
    match = re.match(r"^(\d+)", genotype)
    return match.group(1) if match else genotype


def update_easl_drms(report: list[dict], drm_list: dict, refs: list[dict]) -> list[dict]:
    for report_entry in report:
        sequence_analysis = report_entry.get("sequenceAnalysis")
        if not sequence_analysis:
            continue
        fallback_display = (sequence_analysis.get("strain", {}).get("display", "") or "").replace("Hepatitis C", "").strip()
        genotype = get_genotype_from_subtype(report_entry.get("subtype", ""), fallback_display)
        if not genotype:
            for gene_name in report_entry.get("genes", {}):
                report_entry["genes"][gene_name]["DRMs"] = ""
            continue
        simple_genotype = normalize_genotype(genotype)

        for gene in sequence_analysis.get("alignedGeneSequences") or []:
            gene_name = gene.get("gene", {}).get("name")
            if not gene_name:
                continue
            strict_key = f"{gene_name}_{simple_genotype}"
            strict_gene_drm = drm_list.get(strict_key)
            if not strict_gene_drm:
                if gene_name in report_entry.get("genes", {}):
                    report_entry["genes"][gene_name]["DRMs"] = ""
                continue

            reference_sequence = get_reference_sequence(refs, report_entry.get("subtype", ""), gene_name, fallback_display)
            matched = set()
            for mutation in grouped_pretty_pairwise_mutations(gene):
                position = mutation["position"]
                mut_aa = mutation["mutAA"]
                fallback_ref_aa = mutation["fallbackRefAA"]
                ref_aa = (
                    reference_sequence[position - 1]
                    if reference_sequence and position >= 1 and position <= len(reference_sequence)
                    else fallback_ref_aa
                )
                allowed = strict_gene_drm.get(str(position))
                if not allowed:
                    continue
                if any(char in allowed for char in mut_aa):
                    matched.add(f"{ref_aa}{position}{mut_aa}")
                if mut_aa == "Del" and "#" in allowed:
                    matched.add(f"{ref_aa}{position}Del")

            if gene_name in report_entry.get("genes", {}):
                report_entry["genes"][gene_name]["DRMs"] = ", ".join(
                    sorted(matched, key=lambda text: int(re.search(r"\d+", text).group(0)))
                )
    return report


def analyze(records: list[dict[str, str]], endpoint: str, batch_size: int) -> list[dict]:
    alignments = []
    num_records = len(records)
    progress(f"Aligning: 0 / {num_records} is processed")
    for index in range(0, len(records), batch_size):
        batch = records[index:index + batch_size]
        batch_start = index + 1
        batch_end = index + len(batch)
        progress(f"Aligning: submitting {batch_start}-{batch_end} / {num_records}")
        data = graphql_request(endpoint, ALIGN_QUERY, {"sequences": batch, "geneNames": GENE_NAMES})
        alignments.extend(build_alignment_report(data))
        progress(f"Aligning: {batch_end} / {num_records} is processed")

    validate_gene_exist(alignments)
    validate_no_conflict(alignments)
    alignments = merge_isolates(alignments)
    progress(f"Aligning: merged to {len(alignments)} isolate(s)")

    genotyped = []
    num_alignments = len(alignments)
    progress(f"Genotyping: 0 / {num_alignments} is processed")
    for index in range(0, len(alignments), batch_size):
        batch = alignments[index:index + batch_size]
        batch_start = index + 1
        batch_end = index + len(batch)
        progress(f"Genotyping: submitting {batch_start}-{batch_end} / {num_alignments}")
        genotyped.extend(run_genotyping(batch, endpoint))
        progress(f"Genotyping: {batch_end} / {num_alignments} is processed")

    progress("Calling EASL-extend DRMs")
    drm_list = load_json("genes_drm_easl_extend.json")
    refs = load_json("genes-refs-consensus.json")
    result = update_easl_drms(genotyped, drm_list, refs)
    progress("Calling EASL-extend DRMs: done")
    return result


def build_rows(alignments: list[dict], args: argparse.Namespace) -> list[OrderedDict]:
    rows: list[OrderedDict] = []
    for isolate in alignments:
        genes = isolate["genes"]
        genotype = (isolate.get("genotype", "") or "").upper().replace("GENOTYPE", "").strip()
        subtype = (isolate.get("subtype", "") or "").upper().replace("GENOTYPE", "").strip()

        row = OrderedDict()
        row["IsolateID"] = isolate["IsolateID"]
        row["PersonID"] = f"PID{len(rows) + 1}"
        row["Clone Method"] = args.clone_method
        if args.clone_method != "None":
            row["CloneID"] = ""
        row["NS3 Range"] = genes.get("NS3", {}).get("alignment", {}).get("gene_AA_range", "")
        row["NS5A_NTD Range"] = genes.get("NS5A_NTD", {}).get("alignment", {}).get("gene_AA_range", "")
        row["NS5B Range"] = genes.get("NS5B", {}).get("alignment", {}).get("gene_AA_range", "")
        row["Virus"] = args.organism
        row["Country"] = args.country
        row["Year"] = args.year
        row["Source"] = args.source
        row["Genotype"] = args.genotype or genotype
        row["Subtype"] = args.subtype or subtype
        row["Sequencing Method"] = args.seq_method
        if args.seq_method.lower() != "sanger":
            row["NGS (%)"] = args.ngs_threshold
        row["Treatment"] = args.treatment
        row["DrugClass"] = args.drug_class
        row["PI"] = args.pi
        row["NS5A_I"] = args.ns5a_i
        row["NI"] = args.ni
        row["NNI"] = args.nni
        row["Recent Regimen"] = args.recent_regimen
        row["PI DRMs"] = genes.get("NS3", {}).get("DRMs", "")
        row["NS5A_I DRMs"] = genes.get("NS5A_NTD", {}).get("DRMs", "")
        row["NI DRMs"] = genes.get("NS5B", {}).get("DRMs", "")
        row["NS3 Ins"] = genes.get("NS3", {}).get("insertions", "")
        row["NS5A Ins"] = genes.get("NS5A_NTD", {}).get("insertions", "")
        row["NS5B Ins"] = genes.get("NS5B", {}).get("insertions", "")
        row["Stops"] = "".join(
            f"{gene_name}: {genes[gene_name]['stops']}\n"
            for gene_name in sorted(genes.keys(), key=sort_gene_key)
            if genes[gene_name].get("stops")
        )
        row["DRM Mode"] = args.drm_mode
        rows.append(row)
    return rows


def write_csv_rows(rows: list[OrderedDict], handle) -> None:
    fieldnames = list(rows[0].keys())
    header_writer = csv.writer(handle, quoting=csv.QUOTE_MINIMAL, lineterminator="\n")
    row_writer = csv.writer(handle, quoting=csv.QUOTE_ALL, lineterminator="\n")
    header_writer.writerow(fieldnames)
    for row in rows:
        row_writer.writerow([row.get(fieldname, "") for fieldname in fieldnames])


def write_csv(rows: list[OrderedDict], output: Path | None) -> None:
    if not rows:
        raise ValueError("no rows to write")
    if output is None:
        write_csv_rows(rows, sys.stdout)
        return
    with output.open("w", encoding="utf-8", newline="") as handle:
        write_csv_rows(rows, handle)


def get_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Create the HCV GenBank Step 2 descriptive CSV from a FASTA file."
    )
    parser.add_argument("fasta", type=Path, help="Input FASTA file")
    parser.add_argument("-o", "--output", type=Path, help="Output CSV path. Defaults to stdout.")
    parser.add_argument("--endpoint", default=DEFAULT_ENDPOINT, help="HCV GraphQL endpoint")
    parser.add_argument("--batch-size", type=int, default=20)
    parser.add_argument("--organism", default="HCV")
    parser.add_argument("--source", default="Plasma")
    parser.add_argument("--country", default="USA")
    parser.add_argument("--year", default="", help="Isolate year/date, e.g. 2024 or 2024-05-28")
    parser.add_argument("--clone-method", default="None", choices=["None", "Molecular", "Single genome amplification", "Unique molecular IDs"])
    parser.add_argument("--seq-method", default="Sanger", choices=["Sanger", "Illumina", "Oxford Nanopore", "Ion Torrent", "PacBio"])
    parser.add_argument("--ngs-threshold", default="")
    parser.add_argument("--genotype", default="")
    parser.add_argument("--subtype", default="")
    parser.add_argument("--treatment", default="No", choices=["No", "Yes", "Unknown"])
    parser.add_argument("--drug-class", default="None")
    parser.add_argument("--pi", default="None")
    parser.add_argument("--ns5a-i", default="None")
    parser.add_argument("--ni", default="None")
    parser.add_argument("--nni", default="None")
    parser.add_argument("--recent-regimen", default="None")
    parser.add_argument("--drm-mode", default="EASL-extend")
    return parser


def main() -> int:
    parser = get_parser()
    args = parser.parse_args()
    if args.batch_size < 1:
        parser.exit(1, "error: --batch-size must be >= 1\n")

    try:
        records = parse_fasta(args.fasta)
        alignments = analyze(records, args.endpoint, args.batch_size)
        rows = build_rows(alignments, args)
        write_csv(rows, args.output)
    except (OSError, ValueError, RuntimeError) as exc:
        parser.exit(1, f"error: {exc}\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
