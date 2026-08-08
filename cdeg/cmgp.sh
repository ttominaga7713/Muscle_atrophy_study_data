#!/bin/bash

# =================================================================
# cmgp.sh
# Version: 1.0.0
#
# Comparative Meta-Gene Profiler (CMGP)
# ================================================================= 

VERSION="1.0.0"
IMAGE_NAME="ezojika7713/rnaseq-meta:latest"

# -----------------------------------------------------------------
# Defaults
# -----------------------------------------------------------------
GTF_FILE=""
TARGET_TYPE="all"
FILTER_MODE="include"
PVALUE=0.1
FC_INPUT=""
MODE="edger"
DATASET_CSV=""          
PLOT_TOP_N=5
MAX_MISSING=1.0

# 比較メタプロファイリング用パラメータ
PERMUTATIONS=100000
K_ADJ=0
BOUNDARY_MODE=0
BOUNDARY_VAL="-1.0"
EFPR_MODE=0
EFPR_VAL="-1.0"

SKIP_COUNT=0            
SKIP_EXTRACT=0          
LIST_TYPES=0            

# -----------------------------------------------------------------
# Help and Version
# -----------------------------------------------------------------
show_version() {
    echo "cmgp.sh v${VERSION}"
    exit 0
}

show_help() {
    cat << HELP
Usage: $(basename "$0") [<gtf_file>] [OPTIONS]

Description:
  Comparative Meta-Gene Profiler (CMGP) v${VERSION}.
  Calculates meta-False Discovery Rate (mFDR) using permutations.

Directory Structure & Group Recognition:
  The pipeline automatically scans for '*edgeR_results.tsv' files.
  - Option A (Subdirectories): Place TSV files in subdirectories.
  - Option B (Current Dir): Place all TSV files in the current directory.

Required (unless --rna-types all or --skip-count):
  <gtf_file>               Path to GTF annotation file (.gtf or .gtf.gz)

Gene filtering options:
  -r, --rna-types <str>    Gene type(s) to filter, comma-separated (Default: all)
  -f, --filter <mode>      'include' or 'exclude' (Default: include)
  --pvalue <float>         Raw P-value threshold (T) for single dataset. (Default: 0.1)
  --fc <float>             Fold-change threshold for significance (Default: 1)
  --perms <int>            Number of permutation iterations. (Default: 100000)
  --k-adj <int>            Subtract k from optimal mFDRmin required dataset number (Default: 0)
  --boundary <float>       Find the minimum datasets where mFDR < threshold (e.g., 0.1). Overrides --k-adj.
  --efpr <float>           Find the minimum datasets where eFPR < threshold. Overrides --boundary and --k-adj.

Pipeline control & Extract options:
  -d, --dataset-csv <file> CSV with 'file_name' and 'condition' columns.
  -p, --plot_top <int>     Genes to label in scatter plot (Default: 5)
  --max_missing <float>    Max fraction of missing datasets allowed (Default: 1.0)
  --skip-count             Skip Step 1.
  --skip-extract           Skip Step 2.
  --list-types             Show all gene types and exit.
  -v, --version            Show version information
  -h, --help               Show this help message
HELP
    exit 0
}

# -----------------------------------------------------------------
# Argument parsing
# -----------------------------------------------------------------
if [[ $# -eq 0 ]]; then show_help; fi

while [[ $# -gt 0 ]]; do
    case "$1" in
        -h|--help)          show_help ;;
        -v|--version)       show_version ;;
        -g)                 
                            SPECIES="$2"; VERSION_ARG="$3"
                            if [[ -z "$SPECIES" || -z "$VERSION_ARG" ]]; then
                                echo "Error: -g requires <mouse|human> <version>"
                                exit 1
                            fi
                            SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
                            if [[ -f "$SCRIPT_DIR/count_gene.sh" ]]; then
                                bash "$SCRIPT_DIR/count_gene.sh" -g "$SPECIES" "$VERSION_ARG"
                            else
                                SPECIES_L=$(echo "$SPECIES" | tr '[:upper:]' '[:lower:]')
                                if [[ "$SPECIES_L" == "mouse" ]]; then
                                    URL="https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M${VERSION_ARG}/gencode.vM${VERSION_ARG}.annotation.gtf.gz"
                                    OUT="gencode.vM${VERSION_ARG}.annotation.gtf.gz"
                                elif [[ "$SPECIES_L" == "human" ]]; then
                                    URL="https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_${VERSION_ARG}/gencode.v${VERSION_ARG}.annotation.gtf.gz"
                                    OUT="gencode.v${VERSION_ARG}.annotation.gtf.gz"
                                else
                                    echo "Error: Species must be 'mouse' or 'human'."; exit 1
                                fi
                                echo "Downloading GENCODE GTF..."
                                curl -f -o "$OUT" "$URL" && echo "Downloaded: $OUT" || { echo "Error: Download failed."; rm -f "$OUT"; exit 1; }
                            fi
                            exit 0
                            ;;
        --list-types)       LIST_TYPES=1; shift ;;
        --skip-count)       SKIP_COUNT=1;  shift ;;
        --skip-extract)     SKIP_EXTRACT=1; shift ;;
        -r|--rna-types)     TARGET_TYPE="$2"; shift 2 ;;
        -f|--filter)        FILTER_MODE="$2"; shift 2 ;;
        --pvalue)           PVALUE="$2"; shift 2 ;;
        --fc)               FC_INPUT="$2"; shift 2 ;;
        --perms)            PERMUTATIONS="$2"; shift 2 ;;
        --k-adj)            K_ADJ="$2"; shift 2 ;;
        --boundary)         BOUNDARY_VAL="$2"; BOUNDARY_MODE=1; shift 2 ;;
        --efpr)
            EFPR_MODE=1
            if [[ -n "$2" && "$2" != -* ]]; then
                EFPR_VAL="$2"
                shift 2
            else
                EFPR_VAL="0.05"
                shift 1
            fi
            ;;
        -d|--dataset-csv)   DATASET_CSV="$2"; shift 2 ;;
        -p|--plot_top)      PLOT_TOP_N="$2"; shift 2 ;;
        --venn)             shift ;; 
        --max_missing)      MAX_MISSING="$2"; shift 2 ;;
        -*)
            echo "Error: Unknown option '$1'"
            exit 1
            ;;
        *)
            if [[ -z "$GTF_FILE" ]]; then
                GTF_FILE="$1"
            else
                echo "Error: Unexpected argument '$1'"
                exit 1
            fi
            shift
            ;;
    esac
done

# -----------------------------------------------------------------
# Pre-processing & Validations
# -----------------------------------------------------------------
if [[ $LIST_TYPES -eq 1 ]]; then
    if [[ -z "$GTF_FILE" || ! -f "$GTF_FILE" ]]; then
        echo "Error: GTF file is required to use --list-types."
        echo "Usage: $(basename "$0") <gtf_file> --list-types"
        exit 1
    fi
    echo "================================================================="
    echo " Scanning GTF file for available RNA/Gene types..."
    echo "-----------------------------------------------------------------"
    if [[ "$GTF_FILE" == *.gz ]]; then CAT_CMD="gunzip -c"; else CAT_CMD="cat"; fi
    $CAT_CMD "$GTF_FILE" | awk -F'\t' '
    $3 == "gene" {
        n = split($9, attrs, ";");
        for (i = 1; i <= n; i++) {
            if (index(attrs[i], "gene_type") > 0 || index(attrs[i], "gene_biotype") > 0) {
                split(attrs[i], parts, "\"");
                types[parts[2]] = 1;
            }
        }
    }
    END {
        for (t in types) print "  - " t;
    }' | sort
    echo "================================================================="
    exit 0
fi

if [[ "$K_ADJ" -lt 0 ]]; then
    echo "Error: --k-adj must be >= 0."
    exit 1
fi

if [[ $SKIP_COUNT -eq 0 && "$TARGET_TYPE" != "all" && -z "$GTF_FILE" ]]; then
    echo "Error: GTF file is required when --rna-types is not 'all'."
    exit 1
fi
if [[ $SKIP_COUNT -eq 0 && "$TARGET_TYPE" != "all" && ! -f "$GTF_FILE" ]]; then
    echo "Error: GTF file not found: $GTF_FILE"
    exit 1
fi

if [[ "$MODE" == "edger" ]]; then
    : "${FC_INPUT:=1}"
else
    : "${FC_INPUT:=0}"
fi

if [[ "$TARGET_TYPE" == "all" ]]; then
    FILTER_MODE="all"
fi
SAFE_TARGET=$(echo "$TARGET_TYPE" | tr ',' '_')

if [[ "$TARGET_TYPE" == "all" ]]; then
    OUTPUT_PREFIX="all"
else
    OUTPUT_PREFIX="${SAFE_TARGET}_${FILTER_MODE}"
fi

OUTPUT_DIR="gene_count_${OUTPUT_PREFIX}"
mkdir -p "$OUTPUT_DIR"

echo "================================================================="
echo " Comparative Meta-Gene Profiler (CMGP) v${VERSION}"
echo "-----------------------------------------------------------------"
echo " Target Type     : $TARGET_TYPE ($FILTER_MODE)"
echo " Step 1 (T)      : PValue < $PVALUE, |FC| >= $FC_INPUT"
if [[ $EFPR_MODE -eq 1 ]]; then
    echo " Step 2 (mFDR)   : eFPR Boundary Mode (< $EFPR_VAL Target, ${PERMUTATIONS} Perms)"
elif [[ $BOUNDARY_MODE -eq 1 ]]; then
    echo " Step 2 (mFDR)   : Boundary Mode (< $BOUNDARY_VAL Target, ${PERMUTATIONS} Perms)"
else
    echo " Step 2 (mFDR)   : Baseline Target (< 0.1) (${PERMUTATIONS} Perms, Adj: -$K_ADJ)"
fi
echo "================================================================="

GENE_COUNT_FILES=()

# =================================================================
# STEP 1: Aggregation & Universe Definition
# =================================================================
if [[ $SKIP_COUNT -eq 0 ]]; then
    echo ""
    echo "[ STEP 1 ] Aggregating datasets & Defining Universe per Group"
    
    TEMP_GTF_MAP=$(mktemp /tmp/gtf_map.XXXXXX)
    trap 'rm -f "$TEMP_GTF_MAP"' EXIT

    if [[ "$TARGET_TYPE" == "all" ]]; then
        echo -e "dummy_id\tdummy_type" > "$TEMP_GTF_MAP"
    else
        if [[ "$GTF_FILE" == *.gz ]]; then CAT_CMD="gunzip -c"; else CAT_CMD="cat"; fi
        $CAT_CMD "$GTF_FILE" | awk -F'\t' '
        BEGIN { OFS="\t" }
        $3 == "gene" {
            n = split($9, attrs, ";");
            raw_id = ""; g_type = "";
            for (i = 1; i <= n; i++) {
                split(attrs[i], parts, "\"");
                val = parts[2];
                if (index(attrs[i], "gene_id") > 0) raw_id = val;
                else if (index(attrs[i], "gene_type") > 0 || index(attrs[i], "gene_biotype") > 0) g_type = val;
            }
            if (raw_id != "" && g_type != "") {
                split(raw_id, id_parts, ".");
                clean_id = id_parts[1];
                gsub(/^[ \t]+|[ \t]+$|["\x27\r\n]/, "", clean_id);
                if (clean_id != "") print clean_id "\t" g_type;
            }
        }' > "$TEMP_GTF_MAP"
    fi

    declare -A DIR_LABEL

    CUR_TSVS=( *edgeR_results.tsv )
    if [[ -e "${CUR_TSVS[0]}" ]]; then
        CUR_LABEL=$(basename "$PWD")
        DIR_LABEL["."]="$CUR_LABEL"
    fi

    for subdir in */; do
        subdir="${subdir%/}"
        [[ ! -d "$subdir" ]] && continue
        SUB_TSVS=( "$subdir"/*edgeR_results.tsv )
        if [[ -e "${SUB_TSVS[0]}" ]]; then
            DIR_LABEL["$subdir"]="$subdir"
        fi
    done

    if [[ ${#DIR_LABEL[@]} -eq 0 ]]; then
        echo "  Error: No directories containing *edgeR_results.tsv files found."
        exit 1
    fi

    for dir_path in "${!DIR_LABEL[@]}"; do
        label="${DIR_LABEL[$dir_path]}"
        if [[ "$TARGET_TYPE" == "all" ]]; then
            output_txt="$OUTPUT_DIR/gene_count_all_${label}.txt"
        else
            output_txt="$OUTPUT_DIR/gene_count_${SAFE_TARGET}_${FILTER_MODE}_${label}.txt"
        fi

        if [[ "$dir_path" == "." ]]; then
            tsv_files=( *edgeR_results.tsv )
        else
            tsv_files=( "$dir_path"/*edgeR_results.tsv )
        fi

        file_count=${#tsv_files[@]}
        echo "  Directory : $dir_path ($file_count TSV files)"

        TEMP_EXTRACT=$(mktemp /tmp/extract.XXXXXX)
        TEMP_RESULT=$(mktemp /tmp/result.XXXXXX)
        TEMP_META=$(mktemp /tmp/meta.XXXXXX)

        for f in "${tsv_files[@]}"; do
            [[ ! -f "$f" ]] && continue
            TEMP_EXTRACT_TMP=$(mktemp /tmp/ext_tmp.XXXXXX)

            awk -F'\t' \
                -v p="$PVALUE" -v fc_in="$FC_INPUT" -v mode="$MODE" \
                -v target_str="$TARGET_TYPE" -v fmode="$FILTER_MODE" \
                -v target_all="$( [[ "$TARGET_TYPE" == "all" ]] && echo 1 || echo 0 )" '
            BEGIN {
                id_col=0; sym_col=0; p_col=0; fc_col=0;
                if (!target_all) {
                    split(target_str, t_arr, ",");
                    for (i in t_arr) target_set[t_arr[i]] = 1;
                }
                if (mode == "edger") {
                    if (fc_in < 1) fc_in = 1;
                    fc_cut = log(fc_in) / log(2);
                } else {
                    fc_cut = fc_in;
                }
            }
            NR==FNR { gene_map[$1] = $2; next; }
            FNR==1 {
                for(i=1; i<=NF; i++) {
                    val = $i; sub(/\r$/, "", val); lower = tolower(val);
                    if (mode == "drimseq") {
                        if (lower == "feature_id") id_col=i;
                        else if (lower == "gene_symbol") sym_col=i;
                        else if (lower == "pvalue") p_col=i;
                        else if (lower ~ /^delta_prop/) fc_col=i;
                    } else {
                        if (lower == "geneid") id_col=i;
                        else if (lower == "genesymbol") sym_col=i;
                        else if (lower == "pvalue") p_col=i; 
                        else if (lower == "logfc") fc_col=i;
                    }
                }
            }
            FNR>1 {
                if (p_col > 0 && id_col > 0) {
                    raw_col1 = $id_col;
                    gsub(/^[ \t]+|[ \t]+$|["\x27\r\n]/, "", raw_col1);
                    split(raw_col1, parts, ".");
                    base_id = parts[1];
                    in_map = (base_id in gene_map);
                    is_target = target_all ? 1 : (in_map && gene_map[base_id] in target_set ? 1 : 0);
                    if (target_all || (fmode == "include" && is_target) || (fmode == "exclude" && in_map && !is_target)) {
                        val_p = $p_col;
                        raw_fc = (fc_col > 0) ? $fc_col : 0;
                        if (raw_fc == "") raw_fc = 0;
                        my_key = raw_col1; if (my_key == "") my_key = "-";
                        my_sub = (sym_col > 0) ? $sym_col : "-"; if (my_sub == "") my_sub = "-";
                        is_sig = 0;
                        if (val_p != "" && val_p < p) {
                            abs_fc = (raw_fc < 0) ? -raw_fc : raw_fc;
                            if (abs_fc >= fc_cut) {
                                if (raw_fc > 0) is_sig = 1;
                                else if (raw_fc < 0) is_sig = -1;
                            }
                        }
                        print my_key "\t" my_sub "\t" is_sig "\t" raw_fc;
                    }
                }
            }
            ' "$TEMP_GTF_MAP" "$f" > "$TEMP_EXTRACT_TMP"
            
            cat "$TEMP_EXTRACT_TMP" >> "$TEMP_EXTRACT"
            
            u_c=$(awk -F'\t' '$3==1{c++} END{print c+0}' "$TEMP_EXTRACT_TMP")
            d_c=$(awk -F'\t' '$3==-1{c++} END{print c+0}' "$TEMP_EXTRACT_TMP")
            echo "# Dataset_Hits: UP=$u_c, DOWN=$d_c" >> "$TEMP_META"
            
            rm -f "$TEMP_EXTRACT_TMP"
        done

        awk -F'\t' -v OFS='\t' '
        {
            key = $1; sub_val = $2; status = $3;
            if (!seen[key]) {
                seen[key] = 1; info[key] = (sub_val != "-") ? sub_val : "-";
            } else {
                if (info[key] == "-" && sub_val != "-") info[key] = sub_val;
            }
            if (status == 1)       { up[key]++; }
            else if (status == -1) { down[key]++; }
        }
        END {
            for (k in seen) {
                u = (up[k] != "") ? up[k] : 0;
                d = (down[k] != "") ? down[k] : 0;
                t = u + d;
                print k, info[k], t, u, d;
            }
        }' "$TEMP_EXTRACT" | sort -t$'\t' -k3,3nr > "$TEMP_RESULT"

        if [[ -s "$TEMP_RESULT" ]]; then
            {
                echo "# Dataset_Count: $file_count"
                echo "# Stats Mode: PValue<$PVALUE, FC_Val>=$FC_INPUT"
                cat "$TEMP_META"
                echo -e "Geneid\tGeneSymbol\tTotal_Sig\tSig_UP\tSig_DOWN"
                cat "$TEMP_RESULT"
            } > "$output_txt"
            GENE_COUNT_FILES+=("$output_txt")
        fi
        rm -f "$TEMP_EXTRACT" "$TEMP_RESULT" "$TEMP_META"
    done

else
    echo ""
    echo "[ STEP 1 ] Skipped — collecting existing gene_count_*.txt files"
    if [[ "$TARGET_TYPE" == "all" ]]; then
        COLLECT_PATTERN="gene_count_all_*.txt"
    else
        COLLECT_PATTERN="gene_count_${SAFE_TARGET}_${FILTER_MODE}_*.txt"
    fi
    while IFS= read -r -d '' f; do
        GENE_COUNT_FILES+=("$f")
    done < <(find "$OUTPUT_DIR" . -maxdepth 1 -name "$COLLECT_PATTERN" -print0 2>/dev/null | sort -z)
fi

if [[ $SKIP_EXTRACT -eq 1 ]]; then
    exit 0
fi

# =================================================================
# STEP 2: Meta-Profiling (Permutation & mFDR Calculation)
# =================================================================
echo ""
echo "[ STEP 2 ] Meta-Profiling (Permutation & mFDR Calculation)"
echo "-----------------------------------------------------------------"

INTEGRATED_FILE="$OUTPUT_DIR/.integrated_count_tmp.txt"
TEMP_MAP=$(mktemp /tmp/int_map.XXXXXX)
TEMP_DATA=$(mktemp /tmp/int_data.XXXXXX)
TEMP_COUNTS=$(mktemp /tmp/int_counts.XXXXXX)
TEMP_INT_RESULT=$(mktemp /tmp/int_result.XXXXXX)

for txt_file in "${GENE_COUNT_FILES[@]}"; do
    base=$(basename "$txt_file" .txt)
    if [[ "$TARGET_TYPE" == "all" ]]; then prefix="gene_count_all_"; else prefix="gene_count_${SAFE_TARGET}_${FILTER_MODE}_"; fi
    grp="${base#$prefix}"
    [[ "$grp" == "$base" ]] && grp="$base"
    echo "$grp"$'\t'"$txt_file" >> "$TEMP_MAP"
done

while IFS=$'\t' read -r grp file; do
    awk -F'\t' -v g="$grp" '
    BEGIN { id_col=0; sym_col=0; su_col=0; sd_col=0; ds_count=0; }
    /^# Dataset_Count:/ { split($0, a, ":"); ds_count = a[2]; sub(/^[ \t]+/, "", ds_count); sub(/[ \t\r]+$/, "", ds_count); }
    /^Geneid/ || /^feature_id/ {
        for(i=1; i<=NF; i++) {
            col=$i; sub(/\r$/, "", col);
            if (col == "Geneid" || col == "feature_id") id_col=i;
            if (col == "GeneSymbol" || col == "gene_symbol") sym_col=i;
            if (col == "Sig_UP")   su_col=i;
            if (col == "Sig_DOWN") sd_col=i;
        }
    }
    /^[^#]/ && !/^Geneid/ && !/^feature_id/ {
        if (id_col > 0 && ds_count > 0) {
            my_id  = $id_col;
            my_sym = (sym_col>0) ? $sym_col : "-";
            su = (su_col>0) ? $su_col : 0;
            sd = (sd_col>0) ? $sd_col : 0;
            print my_id "\t" my_sym "\t" g "\t" ds_count "\t" su "\t" sd;
        }
    }
    END { if (ds_count > 0) print g "\t" ds_count >> "'"$TEMP_COUNTS"'" }
    ' "$file" >> "$TEMP_DATA"
done < "$TEMP_MAP"

awk -F'\t' -v OFS='\t' '
NR==FNR {
    group_counts[$1] = $2;
    active_groups[++g_count] = $1;
    next;
}
{
    key = $1; sym = $2; grp = $3; ds_count = $4;
    su = $5; sd = $6;
    if (!seen[key]) { seen[key] = 1; info[key] = sym; }
    sig_up[key, grp]   += su;
    sig_down[key, grp] += sd;
}
END {
    for (k in seen) {
        sum_sig_up=0; sum_sig_down=0;
        hit_details="";
        for (i=1; i<=g_count; i++) {
            g = active_groups[i];
            c_su = sig_up[k,g]+0; c_sd = sig_down[k,g]+0;
            sum_sig_up += c_su; sum_sig_down += c_sd;
            if (c_su>0 || c_sd>0)
                hit_details = hit_details g "(Sig:U" c_su ",D" c_sd ") ";
        }
        sub(/ $/, "", hit_details);
        if (hit_details == "") hit_details = "-";
        
        print k, info[k], sum_sig_up, sum_sig_down, hit_details;
    }
}
' "$TEMP_COUNTS" "$TEMP_DATA" > "$TEMP_INT_RESULT"

{
    echo "# Integrated Count Results"
    echo -n "# Groups Evaluated: "
    awk -F'\t' '{printf "%s(N=%d) ", $1, $2}' "$TEMP_COUNTS"
    echo ""
    echo -e "Geneid\tGeneSymbol\tSig_UP\tSig_DOWN\tGroup_Details"
    cat "$TEMP_INT_RESULT"
} > "$INTEGRATED_FILE"

rm -f "$TEMP_MAP" "$TEMP_DATA" "$TEMP_COUNTS" "$TEMP_INT_RESULT"

if [[ -z "$DATASET_CSV" ]]; then
    DATASET_CSV="$OUTPUT_DIR/dataset_list_auto.csv"
    echo "label,file_name,condition" > "$DATASET_CSV"

    for txt_file in "${GENE_COUNT_FILES[@]}"; do
        base=$(basename "$txt_file" .txt)
        if [[ "$TARGET_TYPE" == "all" ]]; then prefix="gene_count_all_"; else prefix="gene_count_${SAFE_TARGET}_${FILTER_MODE}_"; fi
        grp="${base#$prefix}"
        [[ "$grp" == "$base" ]] && grp="$base"

        if [[ -n "${DIR_LABEL[$grp]+_}" ]]; then dir_path="$grp"; elif [[ -d "$grp" ]]; then dir_path="$grp"; else dir_path="."; fi

        if [[ "$dir_path" == "." ]]; then
            for f in *edgeR_results.tsv; do
                [[ -e "$f" ]] || continue
                fname=$(basename "$f")
                lbl=$(echo "$fname" | grep -oP '^(SRP|ERP|DRP)\d+')
                [[ -z "$lbl" ]] && lbl=$(basename "$f" .tsv)
                echo "$lbl,$f,$grp" >> "$DATASET_CSV"
            done
        else
            for f in "$dir_path"/*edgeR_results.tsv; do
                [[ -e "$f" ]] || continue
                fname=$(basename "$f")
                lbl=$(echo "$fname" | grep -oP '^(SRP|ERP|DRP)\d+')
                [[ -z "$lbl" ]] && lbl=$(basename "$f" .tsv)
                echo "$lbl,$(basename "$f"),$grp" >> "$DATASET_CSV"
            done
        fi
    done
fi

SCRIPT_NAME="temp_core_genes.py"

cat << 'PYEOF' > "$SCRIPT_NAME"
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import seaborn as sns
plt.rcParams['svg.fonttype'] = 'none'
import os
import sys
import hashlib
from pathlib import Path
from functools import reduce
import io
from matplotlib.patches import FancyBboxPatch
import argparse
import re as _re

def simulate_mfdr(univ_size, hits_list, max_votes, perms):
    E_sum = np.zeros(max_votes + 1)
    hits = [h for h in hits_list if h > 0]
    if len(hits) < 1 or univ_size == 0:
        return E_sum
    
    for _ in range(perms):
        votes = np.zeros(univ_size, dtype=np.int32)
        for h in hits:
            idx = np.random.choice(univ_size, h, replace=False)
            votes[idx] += 1
        counts = np.bincount(votes, minlength=max_votes+1)
        for i in range(1, max_votes + 1):
            E_sum[i] += counts[i:].sum()
            
    return E_sum / perms

def get_simulated_E(univ_size, hits_list, max_votes, perms, cache_dir):
    hits = [h for h in hits_list if h > 0]
    if len(hits) < 1 or univ_size == 0:
        return np.zeros(max_votes + 1)
        
    sorted_hits = sorted(hits)
    key_str = f"{univ_size}_{max_votes}_{perms}_{','.join(map(str, sorted_hits))}"
    cache_key = hashlib.md5(key_str.encode('utf-8')).hexdigest()
    
    cache_file = os.path.join(cache_dir, f"sim_{cache_key}.txt")
    if os.path.exists(cache_file):
        print(f"      [Cache] Loaded simulated E values from {cache_key[:8]}...")
        try:
            with open(cache_file, 'r') as f:
                return np.array([float(x) for x in f.read().strip().split(',')])
        except Exception:
            pass
            
    E = simulate_mfdr(univ_size, hits_list, max_votes, perms)
    
    try:
        os.makedirs(cache_dir, exist_ok=True)
        with open(cache_file, 'w') as f:
            f.write(','.join(map(str, E)))
    except Exception:
        pass
        
    return E

def find_k_req(E, N, max_n, univ_size, k_adj, boundary_val, efpr_val):
    mfdrs = {}
    if max_n == 0:
        return -1, 1.0, "N/A"
        
    for i in range(1, max_n + 1):
        mfdrs[i] = (E[i] + 1) / max(N[i], 1)
        
    if efpr_val > 0:
        target_thresh = efpr_val
        for k in sorted(mfdrs.keys()):
            efpr = E[k] / univ_size if univ_size > 0 else 0.0
            if efpr < target_thresh:
                return k, mfdrs[k], "mFDR"
        best_k = max_n
        return best_k, mfdrs[best_k], "mFDR"
        
    elif boundary_val > 0:
        target_thresh = boundary_val
        for k in sorted(mfdrs.keys()):
            m = mfdrs[k]
            if m < target_thresh:
                return k, m, f"mFDR_boundary(<{target_thresh})"
                
        min_m = min(mfdrs.values())
        best_k = min(mfdrs, key=mfdrs.get)
        return best_k, min_m, f"mFDR_boundary(<{target_thresh})"
    else:
        min_m = min(mfdrs.values())
        best_k = min(mfdrs, key=mfdrs.get)
            
        adj_k = max(1, best_k - k_adj)
        adj_m = mfdrs[adj_k]
        
        label = "mFDRmin" if k_adj == 0 else "mFDRmin(adj)"
        
        return adj_k, adj_m, label

def apply_even_tie_breaker(df, base_score_col, logfc_col, max_penalty=1e-9):
    df['sort_logfc'] = np.where(df[base_score_col] >= 0, df[logfc_col], -df[logfc_col])
    df = df.sort_values([base_score_col, 'sort_logfc'], ascending=[False, False])
    
    group_counts = df.groupby(base_score_col).size()
    df['rank_in_group'] = df.groupby(base_score_col).cumcount()
    df['n_in_group'] = df[base_score_col].map(group_counts)
    
    df['penalty'] = (df['rank_in_group'] / np.maximum(1, df['n_in_group'] - 1)) * max_penalty
    is_zero = (df[base_score_col] == 0.0)
    sign = np.sign(df[base_score_col]).replace(0, 1)
    
    final_score = np.where(is_zero,
                           df[logfc_col] * 1e-12,
                           df[base_score_col] - (sign * df['penalty']))
    
    df['GSEA_Score'] = final_score
    return df.drop(columns=['sort_logfc', 'rank_in_group', 'n_in_group', 'penalty'])

def main(integrated_txt, dataset_csv, plot_top_n, max_missing, perms, base_pvalue, fc_input, output_prefix, count_files, base_dir_arg=None, k_adj=0, boundary_val=-1.0, efpr_val=-1.0):
    
    group_thresholds = {}
    cache_dir = ".sim_cache"
    
    print("\n[Simulation] Running Comparative Meta-Gene Profiler (CMGP)...")
    # ------------------------------------------------------------------------
    # STEP 1: Simulate groups & determine mFDR thresholds
    # ------------------------------------------------------------------------
    for cfile in count_files:
        cpath = Path(cfile)
        if not cpath.exists(): continue
        
        basename = cpath.stem
        prefix = f"gene_count_{output_prefix}_"
        grp = basename[len(prefix):] if basename.startswith(prefix) else basename
            
        with open(cpath, 'r') as f:
            lines = f.readlines()
            
        headers = [l for l in lines if l.startswith('#')]
        data_lines = [l for l in lines if not l.startswith('#')]
        
        N_group = 0
        up_hits, down_hits = [], []
        for h in headers:
            if 'Dataset_Count' in h:
                try: N_group = int(h.split(':')[1].strip())
                except: pass
            if 'Dataset_Hits' in h:
                parts = h.split(': ')[1].split(',')
                u = int(parts[0].split('=')[1])
                d = int(parts[1].split('=')[1])
                up_hits.append(u)
                down_hits.append(d)
                
        df_c = pd.read_csv(io.StringIO(''.join(data_lines)), sep='\t')
        univ_size = len(df_c)
        
        N_up, N_down = [0] * (N_group + 1), [0] * (N_group + 1)
        if univ_size > 0:
            for i in range(1, N_group + 1):
                N_up[i] = (df_c['Sig_UP'] >= i).sum()
                N_down[i] = (df_c['Sig_DOWN'] >= i).sum()
                
        print(f"  Group '{grp}' (Universe={univ_size}):")
        E_up = get_simulated_E(univ_size, up_hits, N_group, perms, cache_dir)
        E_down = get_simulated_E(univ_size, down_hits, N_group, perms, cache_dir)
        
        k_up, m_up, lbl_up = find_k_req(E_up, N_up, N_group, univ_size, k_adj, boundary_val, efpr_val)
        k_down, m_down, lbl_down = find_k_req(E_down, N_down, N_group, univ_size, k_adj, boundary_val, efpr_val)

        e_up_val = E_up[k_up] if k_up != -1 else 0.0
        e_down_val = E_down[k_down] if k_down != -1 else 0.0
        efpr_up = e_up_val / univ_size if univ_size > 0 else 0.0
        efpr_down = e_down_val / univ_size if univ_size > 0 else 0.0
        
        group_thresholds[grp] = {
            'n': N_group, 'univ': univ_size,
            'k_up': k_up, 'mFDR_up': m_up, 'lbl_up': lbl_up,
            'E_up_val': e_up_val, 'eFPR_up': efpr_up,
            'k_down': k_down, 'mFDR_down': m_down, 'lbl_down': lbl_down,
            'E_down_val': e_down_val, 'eFPR_down': efpr_down,
            'E_total': e_up_val + e_down_val,
            'eFPR_total': efpr_up + efpr_down
        }
        
        print(f"    UP Req  : {k_up if k_up != -1 else 'N/A'} votes ({lbl_up}={m_up:.4f}) | E={e_up_val:.2f}, eFPR={efpr_up:.4e}")
        print(f"    DOWN Req: {k_down if k_down != -1 else 'N/A'} votes ({lbl_down}={m_down:.4f}) | E={e_down_val:.2f}, eFPR={efpr_down:.4e}")
        
        with open(cpath, 'w') as f:
            for h in headers: f.write(h)
            f.write(f"# --- Comparative Meta-Gene Profiler (CMGP) ({perms} perms) ---\n")
            f.write(f"# Universe Size: {univ_size}\n")
            f.write(f"# UP Required  : {k_up if k_up!=-1 else 'N/A'} ({lbl_up}={m_up:.4f}) | E={e_up_val:.2f}, eFPR={efpr_up:.4e}\n")
            f.write(f"# DOWN Required: {k_down if k_down!=-1 else 'N/A'} ({lbl_down}={m_down:.4f}) | E={e_down_val:.2f}, eFPR={efpr_down:.4e}\n")
            f.write(f"# -----------------------------------------------------\n")
            df_c.to_csv(f, sep='\t', index=False)

    try:
        meta_df = pd.read_csv(dataset_csv)
    except:
        sys.exit(1)
        
    all_csv_files = [os.path.basename(str(row['file_name'])) for _, row in meta_df.iterrows() if not pd.isna(row.get('file_name'))]
    total_datasets = len(all_csv_files)
    
    try: df = pd.read_csv(integrated_txt, sep='\t', comment='#')
    except: sys.exit(1)

    df_integrated = df.copy() 

    condition_map = {}
    for _, row in meta_df.iterrows():
        cond = row['condition'] if 'condition' in meta_df.columns else 'Unknown'
        fname = row['file_name']
        if pd.isna(fname) or pd.isna(cond): continue
        if cond not in condition_map: condition_map[cond] = []
        condition_map[cond].append(os.path.basename(str(fname)))
        
    group_n_map = {cond: len(files) for cond, files in condition_map.items()}
    groups = list(group_n_map.keys())
    n_groups = len(groups)

    # ------------------------------------------------------------------------
    # STEP 2: Evaluate each gene
    # ------------------------------------------------------------------------
    _pat = _re.compile(r'(\S+?)\(Sig:U(\d+),D(\d+)\)')
    
    df_rows = []
    for idx, row in df.iterrows():
        details = str(row.get('Group_Details', ''))
        gene_id = str(row.get('Geneid', ''))
        sym = str(row.get('GeneSymbol', ''))
        if gene_id == "-" or gene_id == "nan" or gene_id == "": continue
        
        recorded = {g: (int(us), int(ds)) for g, us, ds in _pat.findall(details)}
        
        overall_passed = True
        gene_global_dir = None
        ratios = []
        status_recorded = {}
        
        for g in groups:
            n_grp = group_n_map.get(g, 0)
            if n_grp == 0 or g not in group_thresholds:
                overall_passed = False; continue
                
            u_sig, d_sig = recorded.get(g, (0, 0))
            status_recorded[f'{g}_sig_up'] = u_sig
            status_recorded[f'{g}_sig_down'] = d_sig
            
            k_up = group_thresholds[g]['k_up']
            k_down = group_thresholds[g]['k_down']
            
            pass_up = (k_up != -1) and (u_sig >= k_up)
            pass_down = (k_down != -1) and (d_sig >= k_down)
            
            if pass_up and pass_down:
                overall_passed = False
                grp_dir = 'CONFLICT'
                status_recorded[f'{g}_Status'] = 'CONFLICT'
            elif pass_up:
                grp_dir = 'UP'
                status_recorded[f'{g}_Status'] = 'PASS_UP'
            elif pass_down:
                grp_dir = 'DOWN'
                status_recorded[f'{g}_Status'] = 'PASS_DOWN'
            else:
                overall_passed = False
                grp_dir = 'UP' if u_sig > d_sig else 'DOWN'
                status_recorded[f'{g}_Status'] = 'FAIL'
                
            if gene_global_dir is None and grp_dir != 'CONFLICT':
                gene_global_dir = grp_dir
            elif gene_global_dir != grp_dir and overall_passed:
                overall_passed = False 
                
            sig_dom = u_sig if grp_dir == 'UP' else d_sig
            ratios.append(sig_dom / n_grp if n_grp > 0 else 0)
            
        count_ratio = sum(ratios) / n_groups if n_groups > 0 else 0.0
        final_dir = gene_global_dir if overall_passed else ('UP' if row['Sig_UP'] > row['Sig_DOWN'] else 'DOWN')
            
        row_dict = {
            'Geneid': gene_id,
            'GeneSymbol': sym,
            'Dominant_Direction': final_dir,
            'Count_Ratio': count_ratio,
            'is_passed': overall_passed
        }
        row_dict.update(status_recorded)
        df_rows.append(row_dict)
        
    df_eval = pd.DataFrame(df_rows)

    base_dir = Path(base_dir_arg) if base_dir_arg else Path(".")
    
    group_sig_cols = []
    for g in groups:
        group_sig_cols.extend([f'{g}_sig_up', f'{g}_sig_down'])
        
    cols_to_keep = ['Geneid', 'GeneSymbol', 'Dominant_Direction', 'Count_Ratio'] + group_sig_cols + ['is_passed'] + [f'{g}_Status' for g in groups]
    dfs_to_merge = [df_eval[cols_to_keep]]

    all_loaded_files = []
    dataset_stats = []

    for cond, files in condition_map.items():
        for fname in files:
            matched = list(base_dir.rglob(fname))
            if not matched: continue
            file_path = matched[0]
            try:
                tmp = pd.read_csv(file_path, sep='\t')
                gene_col = "GeneSymbol" if "GeneSymbol" in tmp.columns else tmp.columns[0]
                logfc_col = next((col for col in tmp.columns if col.lower() == 'logfc'), None)
                pval_col = next((c for c in tmp.columns if c.lower() in ['pvalue', 'p.value', 'p_value']), None)
                
                sig_count = int(((tmp[pval_col] < base_pvalue) & (tmp[logfc_col].abs() >= fc_input)).sum()) if pval_col and logfc_col else 0
                
                non_samples = {'geneid', 'genesymbol', 'logfc', 'logcpm', 'pvalue', 'fdr', 'feature_id', 'lfcse', 'stat', 'p_value', 'p.value', 'f', 'lr', 'dispersion'}
                sample_cols = [c for c in tmp.columns if c.lower() not in non_samples and not c.lower().startswith('unnamed')]
                n1, n2 = 0, 0
                if len(sample_cols) >= 2:
                    prefixes = [_re.sub(r'[\._-]?\d+$', '', c) for c in sample_cols]
                    from collections import Counter
                    counts = Counter(prefixes).most_common(2)
                    if len(counts) == 2: n1, n2 = counts[0][1], counts[1][1]
                    elif len(counts) == 1 and counts[0][1] >= 2: n1, n2 = counts[0][1] // 2, counts[0][1] - (counts[0][1] // 2)
                    else: n1, n2 = len(sample_cols) // 2, len(sample_cols) - (len(sample_cols) // 2)

                dataset_stats.append({'Dataset': fname.replace('_edgeR_results.tsv', ''), 'Group': cond, 'Sig_Count': sig_count, 'N_Group1': n1, 'N_Group2': n2})

                if logfc_col:
                    sub = tmp[[gene_col, logfc_col]].copy()
                    sub = sub.rename(columns={gene_col: "GeneSymbol", logfc_col: fname})
                    sub = sub.drop_duplicates(subset=["GeneSymbol"]).dropna(subset=["GeneSymbol"])
                    dfs_to_merge.append(sub)
                    all_loaded_files.append(fname)
            except Exception: pass

    merged = reduce(lambda left, right: pd.merge(left, right, on='GeneSymbol', how='left'), dfs_to_merge)

    condition_avg_cols, condition_med_cols = [], []
    for cond, files in condition_map.items():
        valid_files = [f for f in files if f in merged.columns]
        if valid_files:
            cond_avg_col = f"{cond}_Avg_logFC"
            merged[cond_avg_col] = merged[valid_files].mean(axis=1)
            condition_avg_cols.append(cond_avg_col)
            cond_med_col = f"{cond}_Median_logFC"
            merged[cond_med_col] = merged[valid_files].median(axis=1)
            condition_med_cols.append(cond_med_col)

    merged['Final_Avg_logFC'] = merged[condition_avg_cols].mean(axis=1).fillna(0)
    merged['Final_Median_logFC'] = merged[condition_med_cols].median(axis=1)

    # ------------------------------------------------------------------------
    # STEP 3: CORE GENES EXTRACTION AND PLOTTING
    # ------------------------------------------------------------------------
    df_passed = merged[merged['is_passed'] == True].copy()
    
    print(f"\n[Extraction] Target Groups : {n_groups} groups {groups}")
    print(f"  -> Extracted {len(df_passed)} Core Genes that directionally matched in ALL groups.")

    output_txt = f"core_genes_with_logFC_{output_prefix}.txt"

    def write_header():
        ecFPR_up = 1.0
        ecFPR_down = 1.0
        for g in groups:
            ecFPR_up *= group_thresholds[g]['eFPR_up']
            ecFPR_down *= group_thresholds[g]['eFPR_down']
        ecFPR = ecFPR_up + ecFPR_down
        max_univ = max([group_thresholds[g]['univ'] for g in groups]) if groups else 0
        ecFP = ecFPR * max_univ

        with open(output_txt, 'w') as f:
            f.write(f"# --- Comparative Meta-Gene Profiler (CMGP) Results ---\n")
            f.write(f"# Total Datasets (N)   : {total_datasets} across {n_groups} groups {groups}\n")
            f.write(f"# Input Significance   : PValue < {base_pvalue}\n")
            f.write(f"# Permutations         : {perms}\n")
            
            if efpr_val > 0:
                f.write(f"# Mode                 : eFPR Boundary (min datasets where eFPR < {efpr_val})\n")
                f.write(f"# k-adj/boundary Applied: Ignored\n")
            elif boundary_val > 0:
                f.write(f"# Mode                 : Boundary (min datasets where mFDR < {boundary_val})\n")
                f.write(f"# k-adj Applied        : Ignored\n")
            else:
                f.write(f"# Mode                 : Baseline Target (mFDR < 0.1)\n")
                f.write(f"# k-adj Applied        : {k_adj} (subtracted from optimal mFDRmin dataset count)\n")

            f.write(f"# Expected Common FPR (ecFPR)          : {ecFPR:.4e}\n")
            f.write(f"# Base Universe Size (for ecFP)        : {max_univ}\n")
            f.write(f"# Expected Common False Positives (ecFP) : {ecFP:.5f} genes\n")
            f.write(f"# Sorting Criteria     : Count_Ratio (Average of hit ratios across all groups)\n")
            f.write(f"# ----------------------------------------------------\n")
            for g in groups:
                t = group_thresholds[g]
                f.write(f"# Group '{g}' (N={t['n']}, Univ={t['univ']}):\n")
                f.write(f"#    UP Req  : {t['k_up'] if t['k_up']!=-1 else 'N/A'} ({t['lbl_up']}={t['mFDR_up']:.4f}) | E_UP={t['E_up_val']:.2f}, eFPR_UP={t['eFPR_up']:.4e}\n")
                f.write(f"#    DOWN Req: {t['k_down'] if t['k_down']!=-1 else 'N/A'} ({t['lbl_down']}={t['mFDR_down']:.4f}) | E_DOWN={t['E_down_val']:.2f}, eFPR_DOWN={t['eFPR_down']:.4e}\n")
                f.write(f"#    Total Expected Noise (E) : {t['E_total']:.2f} genes, Total eFPR : {t['eFPR_total']:.4e}\n")
            f.write(f"# ----------------------------------------------------\n")

    if df_passed.empty:
        print("  Warning: No genes passed the criteria in all groups. Outputting empty files.")
        write_header()
        empty_cols = ['Geneid', 'GeneSymbol', 'Dominant_Direction', 'Count_Ratio'] + [c for g in groups for c in (f'{g}_sig_up', f'{g}_sig_down')]
        pd.DataFrame(columns=empty_cols).to_csv(output_txt, sep='\t', index=False, mode='a')
        
        plt.figure(figsize=(7, 5))
        plt.text(0.5, 0.5, 'No genes passed the criteria.', ha='center', va='center', fontsize=12)
        plt.axis('off')
        plt.savefig(f"core_genes_scatter_plot_{output_prefix}.png", dpi=300, bbox_inches='tight')
        plt.close()
        return None, group_n_map, df_integrated, group_thresholds

    if max_missing < 1.0:
        file_cols = [f for f in all_loaded_files if f in df_passed.columns]
        missing_counts = (total_datasets - len(file_cols)) + df_passed[file_cols].isna().sum(axis=1)
        missing_frac = missing_counts / total_datasets
        df_passed = df_passed[missing_frac <= max_missing].copy()

    df_passed = df_passed.sort_values(by=['Count_Ratio', 'Final_Avg_logFC'], ascending=[False, False])
    df_passed['Geneid_clean'] = df_passed['Geneid'].apply(lambda x: str(x).split('.')[0])

    score_cols = ['Count_Ratio']
    group_col  = [c for g in groups for c in (f'{g}_sig_up', f'{g}_sig_down')]
    final_cols = ['Final_Avg_logFC', 'Final_Median_logFC']
    cond_stats_cols = condition_avg_cols + condition_med_cols
    status_cols = [f'{g}_Status' for g in groups]
    base_info = ['Geneid', 'Geneid_clean', 'GeneSymbol', 'Dominant_Direction']
    
    placed = set(base_info + score_cols + group_col + final_cols + cond_stats_cols)
    exclude_cols = {'is_passed', 'Group_Details'}.union(status_cols)
    other_cols = [c for c in df_passed.columns if c not in placed and c not in exclude_cols]
    
    df_final = df_passed[base_info + score_cols + group_col + final_cols + cond_stats_cols + other_cols]

    write_header()
    df_final.to_csv(output_txt, sep='\t', index=False, mode='a')

    # ==========================
    # Scatter Plot Drawing
    # ==========================
    plt.figure(figsize=(7, 5))
    ax = plt.gca()
    plot_df = df_final.dropna(subset=['Final_Avg_logFC', 'Count_Ratio']).copy()
    
    if not plot_df.empty:
        sns.scatterplot(
            data=plot_df, x='Final_Avg_logFC', y='Count_Ratio',
            hue='Dominant_Direction', palette={'UP': '#d62728', 'DOWN': '#1f77b4'},
            s=60, alpha=0.7, edgecolor='w', linewidth=0.6, ax=ax
        )
        plt.axhline(0, color='gray', linestyle='-', linewidth=1.2, zorder=1)
        plt.axvline(0, color='gray', linestyle='-', linewidth=1.2, zorder=1)
        ax.xaxis.set_major_locator(ticker.MultipleLocator(2))
        plt.grid(True, linestyle=":", alpha=0.6, zorder=0)
        plt.xlabel("Weighted Average logFC", fontsize=10, fontweight='bold')
        plt.ylabel("Count Ratio", fontsize=10, fontweight='bold')
        
        plt.title("") 
        plt.legend(title="Direction", loc='upper left', bbox_to_anchor=(1.02, 1), fontsize=9)

        x_min, x_max = plot_df['Final_Avg_logFC'].min(), plot_df['Final_Avg_logFC'].max()
        x_range = x_max - x_min
        y_min, y_max = plot_df['Count_Ratio'].min(), plot_df['Count_Ratio'].max()
        y_range = y_max - y_min
        
        plt.xlim(x_min - x_range * 0.35, x_max + x_range * 0.35)
        
        y_min_limit = max(-0.05, y_min - y_range * 0.1)
        ax.set_ylim(y_min_limit, 1.08) 
        
        current_ticks = ax.get_yticks()
        new_ticks = [t for t in current_ticks if t <= 1.0]
        ax.set_yticks(new_ticks)

        if plot_top_n > 0:
            for direction, color in [('UP', 'darkred'), ('DOWN', 'darkblue')]:
                subset = plot_df[plot_df['Dominant_Direction'] == direction]\
                         .sort_values(['Count_Ratio', 'Final_Avg_logFC'], ascending=[False, False if direction=='UP' else True]).head(plot_top_n)
                if subset.empty: continue
                
                subset = subset.sort_values(by='Count_Ratio', ascending=False)
                
                label_start_y = 1.05
                spread_y = max(y_range * 0.35, 0.15)
                label_end_y = max(y_min_limit + 0.1, label_start_y - spread_y)
                
                y_coords = np.linspace(label_start_y, label_end_y, len(subset))
                
                for i, (_, row) in enumerate(subset.iterrows()):
                    tx, ty = row['Final_Avg_logFC'], row['Count_Ratio']
                    offset_x = x_range * 0.18
                    label_x = tx + offset_x if direction == 'UP' else tx - offset_x
                    
                    ax.annotate(row['GeneSymbol'], xy=(tx, ty), xytext=(label_x, y_coords[i]),
                                textcoords='data', ha='left' if direction == 'UP' else 'right',
                                va='center', fontsize=9.5, fontweight='bold', color=color,
                                arrowprops=dict(arrowstyle='-', color='gray', lw=1.2, alpha=0.8))

    plt.savefig(f"core_genes_scatter_plot_{output_prefix}.png", dpi=300, bbox_inches='tight')
    plt.savefig(f"core_genes_scatter_plot_{output_prefix}.svg", format="svg", bbox_inches='tight')
    plt.close()

    return df_final, group_n_map, df_integrated, group_thresholds

def _collect_gene_sets(df_integrated, groups, group_thresholds):
    import re as _re
    _pat = _re.compile(r'(\S+?)\(Sig:U(\d+),D(\d+)\)')
    sets_dict = {g: set() for g in groups}

    for _, row in df_integrated.iterrows():
        gene_id = str(row.get('Geneid', ''))
        if gene_id in ("-", "nan", ""): continue

        details = str(row.get('Group_Details', ''))
        recorded = {g: (int(us), int(ds)) for g, us, ds in _pat.findall(details)}

        pass_up_groups = []
        pass_down_groups = []
        has_conflict = False

        for g in groups:
            if g not in group_thresholds:
                continue
            u_sig, d_sig = recorded.get(g, (0,0))
            k_up = group_thresholds[g]['k_up']
            k_down = group_thresholds[g]['k_down']

            p_up = (k_up != -1) and (u_sig >= k_up)
            p_down = (k_down != -1) and (d_sig >= k_down)

            if p_up and p_down: has_conflict = True
            elif p_up: pass_up_groups.append(g)
            elif p_down: pass_down_groups.append(g)

        if has_conflict:
            continue

        for g in pass_up_groups:
            sets_dict[g].add(f"{gene_id}_UP")
        for g in pass_down_groups:
            sets_dict[g].add(f"{gene_id}_DOWN")

    return sets_dict

def draw_venn(df_integrated, group_n_map, output_prefix, group_thresholds, boundary_val, efpr_val):
    import matplotlib.pyplot as plt
    try:
        from venn import venn
    except ImportError:
        print("  [Venn] Warning: 'venn' library is not installed. Skipping Venn diagram generation.")
        return

    groups = list(group_n_map.keys())
    n_groups = len(groups)
    if not (2 <= n_groups <= 3):
        return

    sets_dict = _collect_gene_sets(df_integrated, groups, group_thresholds)
    labeled_sets = {f"{g}\n(n={group_n_map[g]})": sets_dict[g] for g in groups}

    fig, ax = plt.subplots(figsize=(8, 8))
    venn(labeled_sets, ax=ax)
    
    ax.set_title("Direction-Matched Differential Expression Genes", fontsize=12, fontweight='bold', pad=14)
    plt.savefig(f"venn_diagram_{output_prefix}.png", dpi=300, bbox_inches='tight')
    plt.savefig(f"venn_diagram_{output_prefix}.svg", format="svg", bbox_inches='tight')
    plt.close(fig)

def draw_upset(df_integrated, group_n_map, output_prefix, group_thresholds, boundary_val, efpr_val, target_type, filter_mode):
    try:
        from upsetplot import UpSet, from_contents
    except ImportError:
        print("  [ERROR] 'UpSetPlot' library is missing! RUN pip install UpSetPlot")
        return

    import matplotlib.pyplot as plt
    import matplotlib.ticker as ticker
    import numpy as np
    from matplotlib.patches import FancyBboxPatch

    groups = list(group_n_map.keys())
    sets_dict = _collect_gene_sets(df_integrated, groups, group_thresholds)

    all_genes = set()
    for s in sets_dict.values():
        all_genes |= s

    if not all_genes:
        return

    upset_data = from_contents(sets_dict)

    fig = plt.figure(figsize=(max(8, len(groups) * 1.2), 6))
    upset = UpSet(upset_data, show_counts=True, sort_by='cardinality')
    axes = upset.plot(fig=fig)
    
    axes['intersections'].yaxis.set_major_locator(ticker.MaxNLocator(integer=True))
    axes['intersections'].set_ylabel('Number of Genes', fontsize=10, fontweight='bold')
    
    display_target = "mRNA" if target_type.lower() == "protein_coding" else target_type
    if target_type.lower() == "all":
        title_str = "Model-level DEGs"
    else:
        if filter_mode == "exclude":
            title_str = f"Model-level Non-{display_target} DEGs"
        else:
            title_str = f"Model-level DE{display_target}s"
    plt.suptitle(title_str, fontsize=12, fontweight='bold')
    
    intersections = upset.intersections
    is_all_true = [all(idx) for idx in intersections.index]
    all_true_pos = np.where(is_all_true)[0]

    if len(all_true_pos) > 0:
        pos = all_true_pos[0]
        val = intersections.iloc[pos]
        
        ax_int = axes['intersections']
        ax_mat = axes['matrix']
        
        width_margin = 0.48 
        x_left = pos - width_margin
        x_right = pos + width_margin
        
        ylim_int = ax_int.get_ylim()
        
        y_top_data = val + (ylim_int[1] - ylim_int[0]) * 0.12
        
        p1 = ax_int.transData.transform((x_right, y_top_data))
        p1_fig = fig.transFigure.inverted().transform(p1)
        p0 = ax_int.transData.transform((x_left, 0))
        p0_fig = fig.transFigure.inverted().transform(p0)
        
        mat_bbox = ax_mat.get_position()
        fig_y_bottom = mat_bbox.y0 - 0.015
        
        fig_x = p0_fig[0]
        fig_y = fig_y_bottom
        fig_w = p1_fig[0] - p0_fig[0]
        fig_h = p1_fig[1] - fig_y_bottom
        
        rect = FancyBboxPatch((fig_x, fig_y), fig_w, fig_h,
                              transform=fig.transFigure, boxstyle="round,pad=0.0",
                              fill=False, edgecolor='red', lw=2.0, zorder=100, clip_on=False)
        fig.patches.append(rect)

    plt.savefig(f"upset_plot_{output_prefix}.png", dpi=300, bbox_inches='tight')
    plt.savefig(f"upset_plot_{output_prefix}.svg", format="svg", bbox_inches='tight')
    plt.close('all')

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("integrated_txt")
    parser.add_argument("dataset_csv")
    parser.add_argument("--plot_top",      type=int,   default=5)
    parser.add_argument("--max_missing",   type=float, default=1.0)
    parser.add_argument("--perms",         type=int,   default=100000)
    parser.add_argument("--base_pvalue",   type=float, default=0.1)
    parser.add_argument("--fc_input",      type=float, default=1.0)
    parser.add_argument("--output_prefix", type=str,   default="all")
    parser.add_argument("--do_venn",       action="store_true")
    parser.add_argument("--base_dir",      type=str,   default=None)
    parser.add_argument("--count_files",   nargs="*",  default=[])
    parser.add_argument("--k-adj",         type=int,   default=0)
    parser.add_argument("--boundary_val",  type=float, default=-1.0)
    parser.add_argument("--efpr_val",      type=float, default=-1.0)
    parser.add_argument("--target_type",   type=str,   default="all")
    parser.add_argument("--filter_mode",   type=str,   default="include")
    args = parser.parse_args()
    
    result = main(
        args.integrated_txt, args.dataset_csv, args.plot_top, args.max_missing,
        args.perms, args.base_pvalue, args.fc_input,
        args.output_prefix, args.count_files, base_dir_arg=args.base_dir, k_adj=args.k_adj, boundary_val=args.boundary_val, efpr_val=args.efpr_val
    )
    if args.do_venn and result is not None:
        n_grps = len(result[1])
        if 2 <= n_grps <= 3:
            draw_venn(result[2], result[1], args.output_prefix, result[3], args.boundary_val, args.efpr_val)
        elif n_grps >= 4:
            draw_upset(result[2], result[1], args.output_prefix, result[3], args.boundary_val, args.efpr_val, args.target_type, args.filter_mode)
PYEOF

COUNT_FILES_ARGS=()
for f in "${GENE_COUNT_FILES[@]}"; do
    COUNT_FILES_ARGS+=("/app/$f")
done

BOUNDARY_ARGS=()
if [[ $EFPR_MODE -eq 1 ]]; then
    BOUNDARY_ARGS+=("--efpr_val" "$EFPR_VAL")
elif [[ $BOUNDARY_MODE -eq 1 ]]; then
    BOUNDARY_ARGS+=("--boundary_val" "$BOUNDARY_VAL")
fi

docker run --rm \
    -v "$PWD":/app \
    -w "/app/$OUTPUT_DIR" \
    "$IMAGE_NAME" \
    python3 "/app/$SCRIPT_NAME" "/app/$INTEGRATED_FILE" "/app/$DATASET_CSV" \
        --plot_top "$PLOT_TOP_N" \
        --max_missing "$MAX_MISSING" \
        --perms "$PERMUTATIONS" \
        --base_pvalue "$PVALUE" \
        --fc_input "$FC_INPUT" \
        --output_prefix "$OUTPUT_PREFIX" \
        --base_dir /app \
        --k-adj "$K_ADJ" \
        --do_venn \
        "${BOUNDARY_ARGS[@]}" \
        --target_type "$TARGET_TYPE" \
        --filter_mode "$FILTER_MODE" \
        --count_files "${COUNT_FILES_ARGS[@]}"

rm -f "$SCRIPT_NAME"
rm -f "$INTEGRATED_FILE"

echo ""
echo "================================================================="
echo " Pipeline complete!"
echo " Outputs:"
for f in "${GENE_COUNT_FILES[@]}"; do echo "   [count]   $f"; done
echo "   [core]    $OUTPUT_DIR/core_genes_with_logFC_${OUTPUT_PREFIX}.txt"
echo "   [plot]    $OUTPUT_DIR/core_genes_scatter_plot_${OUTPUT_PREFIX}.png"

if [[ -f "$OUTPUT_DIR/venn_diagram_${OUTPUT_PREFIX}.png" ]]; then
    echo "   [venn]    $OUTPUT_DIR/venn_diagram_${OUTPUT_PREFIX}.png"
fi
if [[ -f "$OUTPUT_DIR/upset_plot_${OUTPUT_PREFIX}.png" ]]; then
    echo "   [upset]   $OUTPUT_DIR/upset_plot_${OUTPUT_PREFIX}.png"
fi

echo "================================================================="